#ifndef CALCMOMENTS_HPP
#define CALCMOMENTS_HPP

#include <chrono>
#include "LBModel.hpp"
#include "Lattice.hpp"
#include <cassert>

class CalcMoments{

private:

    static float _time_calc_moment;

    inline void _get_local_moments(size_t zyx, size_t kdim, const std::vector<int32_t>& c, const float* n, float& rholocal, std::array<float, 3>& ulocal){
        assert(kdim == c.size()/3);
        rholocal = 0.0f;
        ulocal = {0.0f, 0.0f, 0.0f};
        for (auto k=0; k<kdim; ++k){
            auto nk = n[k+zyx*kdim];
            rholocal += nk;
            auto ck = &c[k*3];
            for (auto i=0; i<3; ++i)
                ulocal[i] += nk*ck[i];
        }
        auto rhoinv = 1.0f/rholocal;
        for (auto i=0; i<3; ++i)
            ulocal[i] *= rhoinv;
    }

public:
    CalcMoments(){}

    __attribute__((always_inline))
    inline void operator()(const LBModel *lbmodel, Lattice &lattice){
        auto start = std::chrono::system_clock::now();

        size_t xdim, ydim, zdim, kdim;
        std::tie(zdim, ydim, xdim, kdim) = lattice.n->get_dimensions();

        const auto c = lbmodel->get_lattice_velocities();
        const auto w = lbmodel->get_directional_weights();

        const auto n = lattice.n->get();
        auto rho = lattice.rho->get();
        auto u = lattice.u->get();

        tbb::parallel_for(size_t(1), zdim-1, [this, ydim, xdim, kdim, &c, &w, n, rho, u] (size_t zl){
            for (auto yl=1; yl<ydim-1; ++yl){
                for (auto xl=1; xl<xdim-1; ++xl){
                    float rholocal;
                    std::array<float, 3> ulocal;
                    auto zyx = xl+(yl+zl*ydim)*xdim;
                    _get_local_moments(zyx, kdim, c, n, rholocal, ulocal);
                    rho[zyx] = rholocal;
                    for (auto i=0; i<3; ++i)
                        u[i+zyx*3] = ulocal[i];
                }
            }
        });

        std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
        _time_calc_moment += elapsed.count();
    }

    float get_total_time() const{
        return _time_calc_moment;
    }
};

float CalcMoments::_time_calc_moment = 0.0f;

#endif
