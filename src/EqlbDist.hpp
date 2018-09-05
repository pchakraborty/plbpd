#ifndef SRC_EQLBDIST_HPP_
#define SRC_EQLBDIST_HPP_

#include <array>
#include <cassert>
#include <vector>

#include "LBModel.hpp"

class EqlbDist final {
 private:
    const size_t _kdim;  // number of sub-lattice velocities
    const std::vector<float> &_w;  // directional weights
    const std::vector<std::array<int32_t, 3> > &_c;  // directional velocities

 public:
    explicit EqlbDist(const LBModel *lbmodel)
        : _c(lbmodel->get_directional_velocities()),
          _w(lbmodel->get_directional_weights()),
          _kdim(lbmodel->get_num_directions()) {}
    EqlbDist(const EqlbDist&) = delete;
    EqlbDist& operator=(const EqlbDist&) = delete;
    ~EqlbDist() {}

    __attribute__((always_inline))
    inline void operator()(
        const float rholocal,
        const std::array<float, 3> &ulocal,
        std::vector<float> &nlocal) {
        // start
        auto usq =
            ulocal[0]*ulocal[0] + ulocal[1]*ulocal[1] + ulocal[2]*ulocal[2];
        assert(nlocal.size() == _kdim);
        for (auto k = 0; k < _kdim; ++k) {
            auto cu = _c[k][0]*ulocal[0] + _c[k][1]*ulocal[1] + _c[k][2]*ulocal[2];
            nlocal[k] = _w[k]*rholocal*(1.0+3.0*cu+4.5*cu*cu-1.5*usq);
        }
    }
};

#endif  // SRC_EQLBDIST_HPP_
