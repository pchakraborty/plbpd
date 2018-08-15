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
    : Collision(), _lbmodel(lbmodel), _domain(domain) {
    _tau = 3.0f*_domain->get_fluid_viscosity() + 0.5f;
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
    const auto kdim = simdata.n->get_vector_length();
    const auto c = _lbmodel->get_directional_velocities();
    const auto w = _lbmodel->get_directional_weights();
    const auto extf = _domain->get_external_force();

    std::vector<float> cu(kdim, 0.0f); // scratch space to compute dot(ck, u)

    const auto e = simdata.n->get_extents();
    for (auto zl = e.zbegin; zl < e.zend; ++zl)
        for (auto yl = e.ybegin; yl < e.yend; ++yl)
            for (auto xl = e.xbegin; xl < e.xend; ++xl)
                _collision_kernel(zl, yl, xl, kdim, c, w, extf, simdata, cu.data());
}

// Collision - optimized implementation
void CollisionSRT::_collision_tbb(SimData &simdata) const {
    const auto kdim = simdata.n->get_vector_length();
    const auto c = _lbmodel->get_directional_velocities();
    const auto w = _lbmodel->get_directional_weights();
    const auto f = _domain->get_external_force();

    const auto e = simdata.n->get_extents();
    tbb::parallel_for
    (uint32_t(e.zbegin), e.zend,
     [this, &e, kdim, &c, &w, &f, &simdata] (size_t zl) {
        // cu: scratch space to compute dot(ck,u)
        auto cu = static_cast<float*>(_mm_malloc(kdim*sizeof(float), 64));
        for (auto k = 0; k < kdim; ++k)
            cu[k] = 0.0f;
        for (auto yl = e.ybegin; yl < e.yend; ++yl) {
            for (auto xl = e.xbegin; xl < e.xend; ++xl) {
#if defined(AVX2)
                _collision_kernel_avx2(zl, yl, xl, kdim, c, w, f, simdata, cu);
#else
                _collision_kernel(zl, yl, xl, kdim, c, w, f, simdata, cu);
#endif
            }
        }
        _mm_free(cu);
    });

}

#if defined(AVX2)
// Collision kernel - SIMD (AVX2) implementation
__attribute__((always_inline))
inline void CollisionSRT::_collision_kernel_avx2(
    const size_t zl, const size_t yl, const size_t xl, const size_t kdim,
    const std::vector<int32_t> &c,  // simdata velocities
    const std::vector<float> &w,  // directional weights
    const std::array<float, 3> &ext_force,
    SimData &simdata,
    float * __restrict__ cu) const {  // scratch space to compute dot(ck,u)

    const auto nfpsr = 8;  // (n)umber of (f)loats (p)er (s)imd (r)egister

    // NOTE: I keep getting a segfault if I make these class variables
    const __m256 _one = _mm256_set1_ps(1.0f);
    const __m256 _three = _mm256_set1_ps(3.0f);
    const __m256 _fourPointFive = _mm256_set1_ps(4.5f);
    const __m256 _minusOnePointFive = _mm256_set1_ps(-1.5f);
    const __m256 _omegaVector = _mm256_set1_ps(_omega);
    const __m256 _oneMinusOmega = _mm256_set1_ps(1.0f-_omega);

    // Update u
    auto u_upd = std::array<float, 3>();  // value-initialized to zero
    const float rholocal = simdata.rho->at(zl, yl, xl);
    const float * __restrict__ ulocal = simdata.u->get(zl, yl, xl, 0);
    auto tau_rhoinv = _tau/rholocal;
    for (auto i = 0; i < 3; ++i)
        u_upd[i] = ulocal[i] + ext_force[i]*tau_rhoinv;
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    cu[0] = 0.0f;
    for (auto k = 1; k < kdim; ++k) {
        auto ck = &c[k*3];
        cu[k] = ck[0]*u_upd[0] + ck[1]*u_upd[1] + ck[2]*u_upd[2];
    }

    // nprime(k) = (1.0-_omega)*n(zl,yl,xl,k) + _omega*neq(k);
    // where neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    float * __restrict__ nlocal = simdata.n->get(zl, yl, xl, 0);
    auto _rholocal = _mm256_set1_ps(rholocal);
    auto _kusq = _mm256_set1_ps(-1.5f*usq);  // -1.5*usq
    for (auto k = 0; k < (kdim/nfpsr)*nfpsr; k+=nfpsr) {  // loop unrolling
        auto _w = _mm256_loadu_ps(&w[k]);
        auto _cu = _mm256_load_ps(&cu[k]);
        auto _nk = _mm256_loadu_ps(&nlocal[k]);
        auto _cusq = _mm256_mul_ps(_cu, _cu);
        auto _neq = _mm256_fmadd_ps(_three, _cu, _one);  // _neq = 3*cu(k) + 1
        // _neq += 4.5*cusq(k)
        _neq = _mm256_fmadd_ps(_fourPointFive, _cusq, _neq);
        _neq = _mm256_add_ps(_neq, _kusq);  // _neq += -1.5*usq
        _neq = _mm256_mul_ps(_neq, _rholocal); // _neq *= rholocal
        _neq = _mm256_mul_ps(_neq, _w); // _neq *= w(k)
        // _nk = (1.0-_omega)*n(zl,yl,xl,k) + _omega*_neq
        _nk = _mm256_fmadd_ps(_oneMinusOmega,
                              _nk,
                              _mm256_mul_ps(_omegaVector, _neq));
        _mm256_storeu_ps(nlocal+k, _nk);
    }
    for (auto k = (kdim/nfpsr)*nfpsr; k < kdim; ++k) {  // tail
        auto neq = w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-_omega)*nlocal[k] + _omega*neq;
    }
}
#endif

// Collision kernel - reference implementation
__attribute__((always_inline))
inline void CollisionSRT::_collision_kernel(
    const size_t zl, const size_t yl, const size_t xl, const size_t kdim,
    const std::vector<int32_t> &c,  // directional velocities
    const std::vector<float> &w,  // directional weights
    const std::array<float, 3> &ext_force,
    SimData &simdata,
    float *cu) const {

    // Update u
    auto u_upd = std::array<float, 3>();  // value-initialized to zero
    const float rholocal = simdata.rho->at(zl, yl, xl);
    const float *ulocal = simdata.u->get(zl, yl, xl, 0);
    for (auto i = 0; i < 3; ++i)
        u_upd[i] = ulocal[i] + ext_force[i]*_tau/rholocal;
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    for (auto k = 1; k < kdim; ++k) {
        auto ck = &c[k*3];
        cu[k] = ck[0]*u_upd[0] + ck[1]*u_upd[1] + ck[2]*u_upd[2];
    }

    // neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    // nprime(k) = (1.0-_omega)*n(zl,yl,xl,k) + _omega*neq(k);
    float *nlocal = simdata.n->get(zl, yl, xl, 0);
    for (auto k = 0; k < kdim; ++k) {
        auto neq = w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-_omega)*nlocal[k] + _omega*neq;
    }
}
