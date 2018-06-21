#ifndef CALCMOMENTS_HPP
#define CALCMOMENTS_HPP

#include <chrono>
#include "LBModel.hpp"
#include "SimData.hpp"
#include <cassert>

class CalcMoments final{

private:

    static float _time_calc_moment;

    // A better interface would have been
    // std::pair<float, std::array<float, 3> _get_local_moments(kdim, c, etc)
    // but copying a 3-element std::array at every field node is expensive
    inline void _get_local_moments(
        size_t kdim,
        const std::vector<int32_t>& c,
        const float* nlocal,
        float& rholocal,
        std::array<float, 3>& ulocal) const{
        // Compute density, velocity at a lattice node
        assert(kdim == c.size()/3);
        rholocal = 0.0f;
        ulocal = {0.0f, 0.0f, 0.0f};
        for (auto k=0; k<kdim; ++k){
            auto nk = nlocal[k];
            rholocal += nk;
            auto ck = &c[k*3];
            for (auto i=0; i<3; ++i)
                ulocal[i] += nk*ck[i];
        }
        auto rhoinv = 1.0f/rholocal;
        for (auto it=ulocal.begin(); it!=ulocal.end(); ++it)
            *it *= rhoinv;
    }

public:
    CalcMoments(){}
    CalcMoments(CalcMoments&) = delete;
    CalcMoments& operator=(CalcMoments&) = delete;
    ~CalcMoments(){}
    
    __attribute__((always_inline))
    inline void operator()(const LBModel *lbmodel, SimData &simdata) const{
        auto start = std::chrono::system_clock::now();

        const auto kdim = simdata.n->get_vector_length();
        const auto c = lbmodel->get_directional_velocities();
        const auto w = lbmodel->get_directional_weights();
        const auto e = simdata.n->get_extents();
        
        // NOTE: This implementation (parallelizing the outer loop) is faster
        // than one using tbb::blocked_range3d

        tbb::parallel_for
            (uint32_t(e.zbegin), e.zend, [this, &e, kdim, &simdata, &c, &w] (size_t zl){
                // lambda body - start
                for (auto yl=e.ybegin; yl<e.yend; ++yl){
                    for (auto xl=e.xbegin; xl<e.xend; ++xl){
                        float rholocal;
                        std::array<float, 3> ulocal;
                        auto nlocal = simdata.n->get(zl,yl,xl,0);
                        _get_local_moments(kdim, c, nlocal, rholocal, ulocal);
                        simdata.rho->at(zl,yl,xl) = rholocal;
                        for (auto i=0; i<3; ++i)
                            simdata.u->at(zl,yl,xl,i) = ulocal[i];
                    }
                }
                // lambda body - end
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
