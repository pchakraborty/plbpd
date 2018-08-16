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
      kdim(lbmodel->get_num_directions()),
      c(lbmodel->get_directional_velocities()),
      w(lbmodel->get_directional_weights()),
      extf(domain->get_external_force()) {
    this->tau = 3.0f*domain->get_fluid_viscosity() + 0.5f;
    this->omega = 1./this->tau;
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
    // scratch space to compute dot(ck, u)
    std::vector<float> cu(this->kdim, 0.0f);

    const auto e = simdata.n->get_extents();
    for (auto zl = e.zbegin; zl < e.zend; ++zl)
        for (auto yl = e.ybegin; yl < e.yend; ++yl)
            for (auto xl = e.xbegin; xl < e.xend; ++xl)
                _collision_kernel(zl, yl, xl, simdata, cu.data());
}

// Collision - optimized implementation
void CollisionSRT::_collision_tbb(SimData &simdata) const {
    const auto e = simdata.n->get_extents();

    tbb::parallel_for
    (uint32_t(e.zbegin), e.zend, [this, &e, &simdata] (size_t zl) {
        // thread-local scratch space to compute dot(ck,ueq)
        auto cu = static_cast<float*>(_mm_malloc(this->kdim*sizeof(float), 64));
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

    // // NOTE: Parallelization over z-extent is ~10% faster
    // tbb::parallel_for
    // (tbb::blocked_range3d<uint32_t>
    //  (e.zbegin, e.zend, e.ybegin, e.yend, e.xbegin, e.xend), [this, &simdata]
    //  (const tbb::blocked_range3d<uint32_t> &r) {
    //     // thread-local scratch space to compute dot(ck, ueq)
    //     auto cu = static_cast<float*>(_mm_malloc(this->kdim*sizeof(float), 64));
    //     for (auto zl = r.pages().begin(); zl < r.pages().end(); ++zl)
    //         for (auto yl = r.rows().begin(); yl < r.rows().end(); ++yl)
    //             for (auto xl = r.cols().begin(); xl < r.cols().end(); ++xl) {
    //                 // #if defined(AVX2)
    //                 _collision_kernel_avx2(zl, yl, xl, simdata, cu);
    //                 // #else
    //                 _collision_kernel(zl, yl, xl, simdata, cu);
    //                 // #endif
    //             }
    //     _mm_free(cu);
    // });
}

__attribute__((always_inline))
inline std::array<float, 3> CollisionSRT::_get_updated_u(
    const size_t zl, const size_t yl, const size_t xl,
    const float rholocal, const float *ulocal) const {
    // Update u
    auto u_upd = std::array<float, 3>();  // value-initialized to zero
    auto tau_rhoinv = this->tau/rholocal;
    for (auto i = 0; i < 3; ++i)
        u_upd[i] = ulocal[i] + this->extf[i]*tau_rhoinv;
    return u_upd;
}

__attribute__((always_inline))
inline void CollisionSRT::_get_cu(
    const std::array<float, 3> &u, float* cu) const {
    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    cu[0] = 0.0f;
    for (auto k = 1; k < this->kdim; ++k) {
        auto ck = &this->c[k*3];
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

    const __m256 _one = _mm256_set1_ps(1.0f);
    const __m256 _three = _mm256_set1_ps(3.0f);
    const __m256 _fourPointFive = _mm256_set1_ps(4.5f);
    const __m256 _minusOnePointFive = _mm256_set1_ps(-1.5f);
    const __m256 _omegaVector = _mm256_set1_ps(this->omega);
    const __m256 _oneMinusOmega = _mm256_set1_ps(1.0f-this->omega);

    // Update u and compute usq
    auto u_upd = _get_updated_u(zl, yl, xl, rholocal, ulocal);
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    _get_cu(u_upd, cu);

    // nprime(k) = (1.0-omega)*n(zl,yl,xl,k) + omega*neq(k);
    // where neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    auto _rholocal = _mm256_set1_ps(rholocal);
    auto _kusq = _mm256_set1_ps(-1.5f*usq);  // -1.5*usq
    for (auto k = 0; k < (this->kdim/nfpsr)*nfpsr; k+=nfpsr) {  // unroll loop
        auto _w = _mm256_loadu_ps(&this->w[k]);
        auto _cu = _mm256_load_ps(&cu[k]);
        auto _nk = _mm256_loadu_ps(&nlocal[k]);
        auto _cusq = _mm256_mul_ps(_cu, _cu);
        auto _neq = _mm256_fmadd_ps(_three, _cu, _one);  // _neq = 3*cu(k) + 1
        // _neq += 4.5*cusq(k)
        _neq = _mm256_fmadd_ps(_fourPointFive, _cusq, _neq);
        _neq = _mm256_add_ps(_neq, _kusq);  // _neq += -1.5*usq
        _neq = _mm256_mul_ps(_neq, _rholocal);  // _neq *= rholocal
        _neq = _mm256_mul_ps(_neq, _w);  // _neq *= w(k)
        // _nk = (1.0-omega)*n(zl,yl,xl,k) + omega*_neq
        _nk = _mm256_fmadd_ps(_oneMinusOmega,
                              _nk,
                              _mm256_mul_ps(_omegaVector, _neq));
        _mm256_storeu_ps(nlocal+k, _nk);
    }
    for (auto k = (this->kdim/nfpsr)*nfpsr; k < this->kdim; ++k) {  // tail
        auto neq = this->w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-this->omega)*nlocal[k] + this->omega*neq;
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
    const float * __restrict__ ulocal = simdata.u->get(zl, yl, xl, 0);
    float *nlocal = simdata.n->get(zl, yl, xl, 0);

    // Update u and compute usq
    auto u_upd = _get_updated_u(zl, yl, xl, rholocal, ulocal);
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    _get_cu(u_upd, cu);

    // neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    // nprime(k) = (1.0-omega)*n(zl,yl,xl,k) + omega*neq(k);
    for (auto k = 0; k < this->kdim; ++k) {
        auto neq = this->w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-this->omega)*nlocal[k] + this->omega*neq;
    }
}
