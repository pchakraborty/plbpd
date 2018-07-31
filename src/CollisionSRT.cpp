#include <iostream>
#include <array>
#include <immintrin.h>
#include <mm_malloc.h>
#include <chrono>
#include <algorithm>
#include "tbb/tbb.h"
#include "CollisionSRT.hpp"
#include "EqlbDist.hpp"

CollisionSRT::CollisionSRT(const LBModel *lbmodel, const Domain *domain):
    Collision(), _lbmodel(lbmodel), _domain(domain), _reference(false){
    //_reference = false;
    _tau = 3.0f*_domain->get_fluid_viscosity() + 0.5f;
    _omega = 1./_tau;
}

CollisionSRT::CollisionSRT(const LBModel *lbmodel, const Domain *domain, bool reference):
    Collision(), _lbmodel(lbmodel), _domain(domain), _reference(reference){
    _tau = 3.0f*_domain->get_fluid_viscosity() + 0.5f;
    _omega = 1./_tau;
}

CollisionSRT::~CollisionSRT(){}

void CollisionSRT::operator()(SimData &simdata) const{
    auto start = std::chrono::system_clock::now();

    if (_reference)
        _collision_ref(simdata);
    else
        _collision_tbb(simdata);
    
    std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
    Collision::_time += elapsed.count();
}

// Collision - reference implementation
void CollisionSRT::_collision_ref(SimData &simdata) const{
    const auto kdim = simdata.n->get_vector_length();
    const auto c = _lbmodel->get_directional_velocities();
    const auto w = _lbmodel->get_directional_weights();
    const auto ext_force = _domain->get_external_force();
    const auto e = simdata.n->get_extents();

    for (auto zl=e.zbegin; zl<e.zend; ++zl){
        for (auto yl=e.ybegin; yl<e.yend; ++yl){
            for (auto xl=e.xbegin; xl<e.xend; ++xl){
                auto rholocal = simdata.rho->at(zl,yl,xl);
                auto rhoinv = 1.0f/rholocal;
                std::array<float, 3> ueq;
                for (auto i=0; i<3; ++i)
                    ueq[i] = simdata.u->at(zl,yl,xl,i) + ext_force[i]*_tau*rhoinv;
                auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
                for (auto k=0; k<kdim; ++k){
                    auto ck = &c[k*3];
                    auto cu = ck[0]*ueq[0] + ck[1]*ueq[1] + ck[2]*ueq[2];
                    auto neq = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                    simdata.n->at(zl,yl,xl,k) =
                        (1.0f-_omega)*simdata.n->at(zl,yl,xl,k) + _omega*neq;
                }
            }
        }
    }
}

// Collision - optimized implementation
void CollisionSRT::_collision_tbb(SimData &simdata) const{

    const auto e = simdata.n->get_extents();
    const auto c = _lbmodel->get_directional_velocities();
    const auto w = _lbmodel->get_directional_weights();
    const auto ext_force = _domain->get_external_force();
    const auto kdim = simdata.n->get_vector_length();

    tbb::parallel_for
        (uint32_t(e.zbegin), e.zend,
         [this, &simdata, &e, kdim] (size_t zl){
            const auto c = _lbmodel->get_directional_velocities();
            const auto w = _lbmodel->get_directional_weights();
            const auto ext_force = _domain->get_external_force();
            const float * __restrict__ u = simdata.u->get(0,0,0,0);
            const float * __restrict__ rho = simdata.rho->get(0,0,0);
            float * __restrict__ n = simdata.n->get(0,0,0,0);
            // cu: scratch space to compute dot(ck,u)
            auto cu = static_cast<float*>(_mm_malloc(kdim*sizeof(float), 64));
            for (auto k=0; k<kdim; ++k)
                cu[k] = 0.0f;

            for (auto yl=e.ybegin; yl<e.yend; ++yl){
                for (auto xl=e.xbegin; xl<e.xend; ++xl){
                    auto zyx = simdata.rho->sub2ind(zl,yl,xl);
#if defined(AVX2)
                    _collision_kernel_avx2(zyx, kdim, c, w, ext_force, n, rho, u, cu);
#else
                    _collision_kernel(zyx, kdim, c, w, ext_force, n, rho, u, cu);
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
    const size_t zyx,
    const size_t kdim,
    const std::vector<int32_t> &c, // simdata velocities
    const std::vector<float> &w, // directional weights
    const std::array<float, 3> &ext_force,
    float * __restrict__ n,
    const float * __restrict__ rho,
    const float * __restrict__ u,
    float * __restrict__ cu) const{ // scratch space to compute dot(ck,u)

    const auto fpsr = 8; // number of (f)loats (p)er (s)imd (r)egister

    // Local variables
    auto u_upd = std::array<float, 3>(); // value-initialized to zero
    // Constants
    // NOTE: I keep getting a segfault if I make these class variables
    const __m256 _one = _mm256_set1_ps(1.0f);
    const __m256 _three = _mm256_set1_ps(3.0f);
    const __m256 _fourPointFive = _mm256_set1_ps(4.5f);
    const __m256 _minusOnePointFive = _mm256_set1_ps(-1.5f);
    const __m256 _omegaVector = _mm256_set1_ps(_omega);
    const __m256 _oneMinusOmega = _mm256_set1_ps(1.0f-_omega);

    // Update u
    auto rholocal = rho[zyx];
    auto tau_rhoinv = _tau/rholocal;
    auto ulocal = &u[zyx*3];
    for (auto i=0; i<3; ++i)
        u_upd[i] = ulocal[i] + ext_force[i]*tau_rhoinv;
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    cu[0] = 0.0f;
    for (auto k=1; k<kdim; ++k){
        auto ck = &c[k*3];
        cu[k] = ck[0]*u_upd[0] + ck[1]*u_upd[1] + ck[2]*u_upd[2];
    }

    // nprime(k) = (1.0-_omega)*n(zl,yl,xl,k) + _omega*neq(k);
    // where neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    auto _rholocal = _mm256_set1_ps(rholocal);
    auto _kusq = _mm256_set1_ps(-1.5f*usq); // -1.5*usq
    auto nlocal = &n[zyx*kdim]; // n->get(zl,yl,xl,0)
    for (auto k=0; k<(kdim/fpsr)*fpsr; k+=fpsr){ // loop unrolling
        auto _w = _mm256_loadu_ps(&w[k]);
        auto _cu = _mm256_load_ps(&cu[k]);
        auto _nk = _mm256_loadu_ps(&nlocal[k]);
        auto _cusq = _mm256_mul_ps(_cu, _cu);
        auto _neq = _mm256_fmadd_ps(_three, _cu, _one); // 3*cu(k) + 1
        _neq = _mm256_fmadd_ps(_fourPointFive, _cusq, _neq); // 4.5*cusq(k) + 3*cu(k) + 1
        _neq = _mm256_add_ps(_neq, _kusq); // -1.5*usq + 4.5*cusq(k) + 3*cu(k) + 1
        _neq = _mm256_mul_ps(_neq, _rholocal);
        _neq = _mm256_mul_ps(_neq, _w);
        _nk = _mm256_fmadd_ps(_oneMinusOmega, _nk, _mm256_mul_ps(_omegaVector, _neq));
        _mm256_storeu_ps(nlocal+k, _nk);
    }
    for (auto k=(kdim/fpsr)*fpsr; k<kdim; ++k){ // tail
        auto neq = w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-_omega)*nlocal[k] + _omega*neq;
    }
}
#endif

// Collision kernel - reference implementation
__attribute__((always_inline))
inline void CollisionSRT::_collision_kernel(
    const size_t zyx,
    const size_t kdim,
    const std::vector<int32_t> &c, // directional velocities
    const std::vector<float> &w, // directional weights
    const std::array<float, 3> &ext_force,
    float * __restrict__ n,
    const float * __restrict__ rho,
    const float * __restrict__ u,
    float *cu) const{

    // Local variables
    auto u_upd = std::array<float, 3>(); // value-initialized to zero

    // Update u
    auto rholocal = rho[zyx];
    auto tau_rhoinv = _tau/rholocal;
    for (auto i=0; i<3; ++i)
        u_upd[i] = u[zyx*3+i] + ext_force[i]*tau_rhoinv;
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    for (auto k=1; k<kdim; ++k){
        auto ck = &c[k*3];
        cu[k] = ck[0]*u_upd[0] + ck[1]*u_upd[1] + ck[2]*u_upd[2];
    }

    // neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    // nprime(k) = (1.0-_omega)*n(zl,yl,xl,k) + _omega*neq(k);
    auto nlocal = &n[zyx*kdim+0]; // n->get(zl,yl,xl,0)
    for (auto k=0; k<kdim; ++k){
        auto neq = w[k]*rholocal*(1.0+3.0*cu[k]+4.5*cu[k]*cu[k]-1.5*usq);
        nlocal[k] = (1.0f-_omega)*nlocal[k] + _omega*neq;
    }
}