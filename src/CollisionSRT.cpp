#include <immintrin.h>
#include <mm_malloc.h>
#include <chrono>
#include <array>
#include <vector>
#include <algorithm>
#include "tbb/tbb.h"

#include "CollisionSRT.hpp"
#include "EqlbDist.hpp"

namespace chrono = std::chrono;

CollisionSRT::CollisionSRT(const LBModel *lbmodel, const Domain *domain)
    : Collision(),
      _c(lbmodel->get_directional_velocities()),
      _w(lbmodel->get_directional_weights()),
      _extf(domain->get_external_force()) {
    _tau = 3.0f*domain->get_fluid_viscosity() + 0.5f;
    _omega = 1./_tau;
}

CollisionSRT::~CollisionSRT() {}

void CollisionSRT::operator()(SimData &simdata, bool reference) const {
    auto start = chrono::system_clock::now();
    if (reference)
        _collision_ref(simdata);
    else
        _collision_tbb(simdata);
    chrono::duration<float> elapsed = chrono::system_clock::now()-start;
    Collision::_time += elapsed.count();
}

// Collision - reference implementation
void CollisionSRT::_collision_ref(SimData &simdata) const {
    auto kdim = simdata.n->get_vector_length();
    // scratch space to compute dot(ck, u)
    std::vector<float> cu(kdim, 0.0f);

    const auto e = simdata.n->get_extents();
    for (auto zl = e.zbegin; zl < e.zend; ++zl)
        for (auto yl = e.ybegin; yl < e.yend; ++yl)
            for (auto xl = e.xbegin; xl < e.xend; ++xl)
                _collision_kernel(zl, yl, xl, simdata, cu.data());
}

// Collision - optimized implementation
void CollisionSRT::_collision_tbb(SimData &simdata) const {
    const auto e = simdata.n->get_extents();
    // NOTE: In this particular case, tbb::parallel_for over
    // z-extent is ~10% faster than that over tbb::blocked_range3d
    tbb::parallel_for
    (uint32_t(e.zbegin), e.zend, [this, &e, &simdata] (size_t zl) {
        auto kdim = simdata.n->get_vector_length();
        // thread-local scratch space to compute dot(ck,ueq)
        auto cu = static_cast<float*>(_mm_malloc(kdim*sizeof(float), 64));
        for (auto yl = e.ybegin; yl < e.yend; ++yl) {
            for (auto xl = e.xbegin; xl < e.xend; ++xl) {
#if defined(AVX2)
                _collision_kernel_avx2(zl, yl, xl, simdata, cu);
#else
                _collision_kernel(zl, yl, xl, simdata, cu);
#endif
            }
        }
        _mm_free(cu);
    });
}

__attribute__((always_inline))
inline std::array<float, 3> CollisionSRT::_get_updated_u(
    const size_t zl, const size_t yl, const size_t xl,
    const float rholocal, const float *ulocal) const {
    // Update u
    auto u_upd = std::array<float, 3>();  // value-initialized to zero
    auto tau_rhoinv = _tau/rholocal;
    for (auto i = 0; i < 3; ++i)
        u_upd[i] = ulocal[i] + _extf[i]*tau_rhoinv;
    return u_upd;
}

__attribute__((always_inline))
inline void CollisionSRT::_get_cu(
    const std::array<float, 3> &u,
    float *cu, const size_t cu_len) const {
    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    cu[0] = 0.0f;
    for (auto k = 1; k < cu_len; ++k) {
        auto ck = &_c[k*3];
        cu[k] = ck[0]*u[0] + ck[1]*u[1] + ck[2]*u[2];
    }
}

#if defined(AVX2)
// Collision kernel - SIMD (AVX2) implementation
__attribute__((always_inline))
inline void CollisionSRT::_collision_kernel_avx2(
    const size_t zl, const size_t yl, const size_t xl,
    SimData &simdata,
    float * __restrict__ cu) const {  // scratch space to compute dot(ck,u)

    const auto nfpsr = 8;  // (n)umber of (f)loats (p)er (s)imd (r)egister

    const float rholocal = simdata.rho->at(zl, yl, xl);
    const float * __restrict__ ulocal = simdata.u->get(zl, yl, xl, 0);
    float * __restrict__ nlocal = simdata.n->get(zl, yl, xl, 0);
    const auto kdim = simdata.n->get_vector_length();

    const __m256 m_one = _mm256_set1_ps(1.0f);
    const __m256 m_three = _mm256_set1_ps(3.0f);
    const __m256 m_four_point_five = _mm256_set1_ps(4.5f);
    const __m256 m_minus_one_point_five = _mm256_set1_ps(-1.5f);
    const __m256 m_omega = _mm256_set1_ps(_omega);
    const __m256 m_one_minus_omega = _mm256_set1_ps(1.0f-_omega);

    // Update u and compute usq
    auto u_upd = _get_updated_u(zl, yl, xl, rholocal, ulocal);
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    _get_cu(u_upd, cu, kdim);

    // nprime(k) = (1.0-omega)*n(zl,yl,xl,k) + omega*neq(k);
    // where neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    auto m_rholocal = _mm256_set1_ps(rholocal);
    auto m_kusq = _mm256_set1_ps(-1.5f*usq);  // -1.5*usq
    for (auto k = 0; k < (kdim/nfpsr)*nfpsr; k+=nfpsr) {  // unroll loop
        auto m_w = _mm256_loadu_ps(&_w[k]);
        auto m_cu = _mm256_load_ps(&cu[k]);
        auto m_nk = _mm256_loadu_ps(&nlocal[k]);
        auto m_cusq = _mm256_mul_ps(m_cu, m_cu);
        // neq = 3*cu(k) + 1
        auto m_neq = _mm256_fmadd_ps(m_three, m_cu, m_one);
        // neq += 4.5*cusq(k)
        m_neq = _mm256_fmadd_ps(m_four_point_five, m_cusq, m_neq);
        m_neq = _mm256_add_ps(m_neq, m_kusq);  // neq += -1.5*usq
        m_neq = _mm256_mul_ps(m_neq, m_rholocal);  // neq *= rholocal
        m_neq = _mm256_mul_ps(m_neq, m_w);  // neq *= w(k)
        // nk = (1.0-omega)*n(zl,yl,xl,k) + omega*neq
        m_nk = _mm256_fmadd_ps(m_one_minus_omega,
                               m_nk,
                               _mm256_mul_ps(m_omega, m_neq));
        _mm256_storeu_ps(nlocal+k, m_nk);
    }
    for (auto k = (kdim/nfpsr)*nfpsr; k < kdim; ++k) {  // tail
        auto neq = _w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-_omega)*nlocal[k] + _omega*neq;
    }
}
#endif

// Collision kernel - reference implementation
__attribute__((always_inline))
inline void CollisionSRT::_collision_kernel(
    const size_t zl, const size_t yl, const size_t xl,
    SimData &simdata,
    float *cu) const {

    const float rholocal = simdata.rho->at(zl, yl, xl);
    const float *ulocal = simdata.u->get(zl, yl, xl, 0);
    float *nlocal = simdata.n->get(zl, yl, xl, 0);
    const auto kdim = simdata.n->get_vector_length();

    // Update u and compute usq
    auto u_upd = _get_updated_u(zl, yl, xl, rholocal, ulocal);
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    _get_cu(u_upd, cu, kdim);

    // neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    // nprime(k) = (1.0-omega)*n(zl,yl,xl,k) + omega*neq(k);
    for (auto k = 0; k < kdim; ++k) {
        auto neq = _w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-_omega)*nlocal[k] + _omega*neq;
    }
}
