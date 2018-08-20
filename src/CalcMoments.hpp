#ifndef CALCMOMENTS_HPP
#define CALCMOMENTS_HPP

#include <chrono>
#include <vector>
#include <array>
#include "tbb/tbb.h"

#include "LBModel.hpp"
#include "SimData.hpp"

namespace chrono = std::chrono;

class CalcMoments final {
 private:
    static float _time_calc_moment;
    const std::vector<std::array<int32_t, 3> > &_c;
    const std::vector<float> &_w;
    const size_t _kdim;

    // A better interface would have been
    // std::pair<float, std::array<float, 3> _get_local_moments(kdim, c, etc)
    // but copying a 3-element std::array at every field node is expensive
    inline void _get_local_moments(
        const float* nlocal,
        float& rholocal,
        std::array<float, 3>& ulocal) const {
        // Compute density, velocity at a lattice node
        assert(_kdim == _c.size()/3);
        rholocal = 0.0f;
        ulocal = {0.0f, 0.0f, 0.0f};
        for (auto k = 0; k < _kdim; ++k) {
            auto nk = nlocal[k];
            rholocal += nk;
            for (auto i = 0; i < 3; ++i)
                ulocal[i] += nk*_c[k][i];
        }
        auto rhoinv = 1.0f/rholocal;
        for (auto it = ulocal.begin(); it != ulocal.end(); ++it)
            *it *= rhoinv;
    }

 public:
    explicit CalcMoments(const LBModel *lbmodel);
    CalcMoments(CalcMoments&) = delete;
    CalcMoments& operator=(CalcMoments&) = delete;
    ~CalcMoments();
    float get_total_time() const;

    __attribute__((always_inline))
    inline void operator()(SimData &simdata) const {
        auto start = chrono::system_clock::now();

        // NOTE: This implementation (parallelizing the outer loop) is faster
        // than one using tbb::blocked_range3d

        const auto e = simdata.n->get_extents();
        tbb::parallel_for
        (uint32_t(e.zbegin), e.zend, [this, &e, &simdata] (size_t zl) {
            // lambda body - start
            for (auto yl = e.ybegin; yl < e.yend; ++yl) {
                for (auto xl = e.xbegin; xl < e.xend; ++xl) {
                    float rholocal;
                    std::array<float, 3> ulocal;
                    auto nlocal = simdata.n->get(zl, yl, xl, 0);
                    _get_local_moments(nlocal, rholocal, ulocal);
                    simdata.rho->at(zl, yl, xl) = rholocal;
                    for (auto i = 0; i < 3; ++i)
                        simdata.u->at(zl, yl, xl, i) = ulocal[i];
                }
            }
            // lambda body - end
        });

        chrono::duration<float> elapsed = chrono::system_clock::now()-start;
        _time_calc_moment += elapsed.count();
    }
};

#endif
