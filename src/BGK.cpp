#include <iostream>
#include "BGK.hpp"
#include "tbb/tbb.h"
#include <array>
#include <immintrin.h>
#include <mm_malloc.h>
#include "EqlbDist.hpp"
#include <chrono>

BGK::BGK(const LBModel *lbmodel, const Domain *domain):
    LBDynamics(), _lbmodel(lbmodel), _domain(domain) {
    _tau = 3.0f*_domain->getFluidViscosity() + 0.5f;
    _omega = 1./_tau;
}

BGK::~BGK(){}

void BGK::collideAndStream(Lattice &lattice) const{
    auto start = std::chrono::system_clock::now();

    // // Reference implementations
    // _collide_ref(lattice);
    // _stream_ref(lattice);

    // Optimized implementations
    _collide(lattice);
    _stream(lattice);

    std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
    LBDynamics::_timeTakenByCollideAndStream += elapsed.count();
}

void BGK::_stream_ref(Lattice &lattice) const{
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto c = _lbmodel->getLatticeVelocities();
    array4f __restrict__ *n = lattice.n;
    array4f __restrict__ *ntmp = lattice.ntmp;
    for(auto zl=1; zl<zdim+1; ++zl){
        for (auto yl=1; yl<ydim+1; ++yl){
            for (auto xl=1; xl<xdim+1; ++xl){
                for (auto k=0; k<kdim; ++k){
                    auto ck = &c[k*3];
                    ntmp->at(zl+ck[2],yl+ck[1],xl+ck[0],k) = n->at(zl,yl,xl,k);
                }
            }
        }
    }
    // swap n and ntmp
    array4f *tmp = lattice.n;
    lattice.n = lattice.ntmp;
    lattice.ntmp = tmp;
}

void BGK::_stream(Lattice &lattice) const{
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    tbb::parallel_for(size_t(1), zdim+1, [this, &lattice, ydim, xdim] (size_t zl){
        const auto kdim = _lbmodel->getNumberOfDirections();
        const auto c = _lbmodel->getLatticeVelocities();
        array4f __restrict__ *n = lattice.n;
        array4f __restrict__ *ntmp = lattice.ntmp;
        for (auto yl=1; yl<ydim+1; ++yl){
            for (auto xl=1; xl<xdim+1; ++xl){
                for (auto k=0; k<kdim; ++k){
                    auto ck = &c[k*3];
                    ntmp->at(zl+ck[2],yl+ck[1],xl+ck[0],k) = n->at(zl,yl,xl,k);
                }
            }
        }
    });
    // swap n and ntmp
    array4f *tmp = lattice.n;
    lattice.n = lattice.ntmp;
    lattice.ntmp = tmp;
}

// Collision - reference implementation
void BGK::_collide_ref(Lattice &lattice) const{
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();
    const auto extForce = _domain->getExternalForce();
    for (auto zl=1; zl<zdim+1; ++zl){
        for (auto yl=1; yl<ydim+1; ++yl){
            for (auto xl=1; xl<xdim+1; ++xl){
                auto rholocal = lattice.rho->at(zl,yl,xl);
                auto rhoinv = 1.0f/rholocal;
                std::array<float, 3> ueq;
                for (auto i=0; i<3; ++i)
                    ueq[i] = lattice.u->at(zl,yl,xl,i) + extForce[i]*_tau*rhoinv;
                auto usq = ueq[0]*ueq[0] + ueq[1]*ueq[1] + ueq[2]*ueq[2];
                for (auto k=0; k<kdim; ++k){
                    auto ck = &c[k*3];
                    auto cu = ck[0]*ueq[0] + ck[1]*ueq[1] + ck[2]*ueq[2];
                    auto neq = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
                    lattice.n->at(zl,yl,xl,k) =
                        (1.0f-_omega)*lattice.n->at(zl,yl,xl,k) + _omega*neq;
                }
            }
        }
    }
}

// Collision - optimized implementation
void BGK::_collide(Lattice &lattice) const{

    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();

    //for (auto zl=1; zl<zdim+1; ++zl){
    tbb::parallel_for(size_t(1), zdim+1, [this, &lattice, ydim, xdim] (size_t zl){
        const auto kdim = _lbmodel->getNumberOfDirections();
        const auto c = _lbmodel->getLatticeVelocities();
        const auto w = _lbmodel->getDirectionalWeights();
        const auto extForce = _domain->getExternalForce();
        // cu is a scratch vector to compute dot(ck,u)
        auto cu = static_cast<float*>(_mm_malloc(kdim*sizeof(float), 64));
        for (auto k=0; k<kdim; ++k) cu[k] = 0.0f;
        float * __restrict__ u = lattice.u->get();
        float * __restrict__ rho = lattice.rho->get();
        float * __restrict__ n = lattice.n->get();
        for (auto yl=1; yl<ydim+1; ++yl){
            for (auto xl=1; xl<xdim+1; ++xl){
                auto ndx3d = xl+(yl+zl*(ydim+2))*(xdim+2);
                _collide_kernel(ndx3d, kdim, c, w, extForce, n, rho, u, cu);
            }
        }
        _mm_free(cu);
    });

}

// Collision kernel - SIMD (AVX2) implementation
__attribute__((always_inline))
inline void BGK::_collide_kernel(
    const size_t zyx,
    const size_t kdim,
    const std::vector<int32_t> &c, // lattice velocities
    const std::vector<float> &w, // directional weights
    const std::array<float, 3> &extForce,
    float * __restrict__ n,
    const float * __restrict__ rho,
    const float * __restrict__ u,
    float *cu) const{

    const auto fpsr = 8; // number of (f)loats (p)er (s)imd (r)egister

    // Local variables
    auto u_upd = std::array<float, 3>(); // value-initialized to zero
    // Constants - I keep getting a segfault if I make these class variables
    const __m256 _one = _mm256_set1_ps(1.0f);
    const __m256 _three = _mm256_set1_ps(3.0f);
    const __m256 _fourPointFive = _mm256_set1_ps(4.5f);
    const __m256 _minusOnePointFive = _mm256_set1_ps(-1.5f);
    const __m256 _omegaVector = _mm256_set1_ps(_omega);
    const __m256 _oneMinusOmega = _mm256_set1_ps(1.0f-_omega);

    // Update u
    auto rholocal = rho[zyx];
    auto tau_rhoinv = _tau/rholocal;
    for (auto i=0; i<3; ++i)
        u_upd[i] = u[zyx*3+i] + extForce[i]*tau_rhoinv;
    auto usq = u_upd[0]*u_upd[0] + u_upd[1]*u_upd[1] + u_upd[2]*u_upd[2];

    // cu(k) = c(k,i)*ueq(i), 19 3x3 dot products for D3Q19
    for (auto k=1; k<kdim; ++k){
        auto ck = &c[k*3];
        cu[k] = ck[0]*u_upd[0] + ck[1]*u_upd[1] + ck[2]*u_upd[2];
    }

    // neq(k) = w(k)*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
    // nprime(k) = (1.0-_omega)*n(zl,yl,xl,k) + _omega*neq(k);
    auto _rholocal = _mm256_set1_ps(rholocal);
    auto _kusq = _mm256_set1_ps(-1.5f*usq); // -1.5*usq
    auto nlocal = &n[zyx*kdim+0]; // n->get(zl,yl,xl,0)
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

void BGK::calcMoments(Lattice &lattice) const{
    auto start = std::chrono::system_clock::now();

    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = _domain->getDimensions();
    const auto kdim = _lbmodel->getNumberOfDirections();
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();

    // const array4f *n = lattice.n;
    const auto n = lattice.n->get();
    auto rho = lattice.rho->get();
    auto u = lattice.u->get();

    tbb::parallel_for(size_t(1), zdim+1, [this, ydim, xdim, kdim, &c, &w, n, rho, u] (size_t zl){
        for (auto yl=1; yl<ydim+1; ++yl){
            for (auto xl=1; xl<xdim+1; ++xl){
                auto ndx3d = xl+(yl+zl*(ydim+2))*(xdim+2);
                float rholocal = 0.0f;
                std::array<float, 3> ulocal = {0.0f, 0.0f, 0.0f};
                for (auto k=0; k<kdim; ++k){
                    auto nk = n[k+ndx3d*kdim];
                    rholocal += nk; // \sum_k n[k]
                    auto ck = &c[k*3];
                    for (auto i=0; i<3; ++i)
                        ulocal[i] += nk*ck[i]; // \sum_k n[k]*c[k][i], i=0,1,2
                }
                rho[ndx3d] = rholocal;
                auto rhoinv = 1.0f/rholocal;
                for (auto i=0; i<3; ++i)
                    u[i+ndx3d*3] = ulocal[i]*rhoinv;
            }
        }
    });

    std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
    LBDynamics::_timeTakenByCalcMoments += elapsed.count();
}
