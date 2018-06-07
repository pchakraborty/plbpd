#ifndef CALCMOMENTS_HPP
#define CALCMOMENTS_HPP

#include <chrono>
#include "LBModel.hpp"
#include "Lattice.hpp"

class CalcMoments{

private:
    static float _time_calc_moment;

public:
    CalcMoments(){}

    __attribute__((always_inline))
    inline void operator()(const LBModel *lbmodel, Lattice &lattice){
        auto start = std::chrono::system_clock::now();

        size_t xdim, ydim, zdim, kdim;
        std::tie(zdim, ydim, xdim, kdim) = lattice.n->getDimensions();

        const auto c = lbmodel->getLatticeVelocities();
        const auto w = lbmodel->getDirectionalWeights();

        const auto n = lattice.n->get();
        auto rho = lattice.rho->get();
        auto u = lattice.u->get();
        
        tbb::parallel_for(size_t(1), zdim-1, [ydim, xdim, kdim, &c, &w, n, rho, u] (size_t zl){
            for (auto yl=1; yl<ydim-1; ++yl){
                for (auto xl=1; xl<xdim-1; ++xl){
                    auto ndx3d = xl+(yl+zl*ydim)*xdim;
                    auto rholocal = 0.0f;
                    std::array<float, 3> ulocal = {0.0f, 0.0f, 0.0f};
                    for (auto k=0; k<kdim; ++k){
                        auto nk = n[k+ndx3d*kdim];
                        rholocal += nk;
                        auto ck = &c[k*3];
                        for (auto i=0; i<3; ++i)
                            ulocal[i] += nk*ck[i];
                    }
                    rho[ndx3d] = rholocal;
                    auto rhoinv = 1.0f/rholocal;
                    for (auto i=0; i<3; ++i)
                        u[i+ndx3d*3] = ulocal[i]*rhoinv;
                }
            }
        });

        std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
        _time_calc_moment += elapsed.count();
    }

    float get_time_taken() const{
        return _time_calc_moment;
    }
};

float CalcMoments::_time_calc_moment = 0.0f;

#endif
