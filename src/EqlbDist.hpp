#ifndef EQLBDIST_HPP
#define EQLBDIST_HPP

#include <array>
#include <cassert>
#include <vector>

#include "LBModel.hpp"

class EqlbDist final {
 public:
    EqlbDist() {}
    // EqlbDist(EqlbDist&) = delete;
    EqlbDist& operator=(EqlbDist&) = delete;
    ~EqlbDist() {}

    __attribute__((always_inline))
    inline void operator()(
        const LBModel *lbmodel,
        const float rholocal,
        const std::array<float, 3> &ulocal,
        std::vector<float> &nlocal) {
        // start
        const auto kdim = lbmodel->get_num_directions();
        const auto w = lbmodel->get_directional_weights();
        const auto c = lbmodel->get_directional_velocities();
        auto usq =
            ulocal[0]*ulocal[0] + ulocal[1]*ulocal[1] + ulocal[2]*ulocal[2];
        assert(nlocal.size() == kdim);
        for (auto k = 0; k < kdim; ++k) {
            auto ck = &c[k*3];
            auto cu = ck[0]*ulocal[0] + ck[1]*ulocal[1] + ck[2]*ulocal[2];
            nlocal[k] = w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
        }
    }
};

#endif
