#ifndef COLLISION_SRT_HPP
#define COLLISION_SRT_HPP

#include <vector>
#include <array>

#include "Collision.hpp"
#include "LBModel.hpp"
#include "Domain.hpp"

class CollisionSRT final: public Collision {
 private:
    // It's convenient to store the input const references
    const size_t kdim;  // number of sub-lattice directions
    const std::vector<int32_t> &c;  // directional velocities
    const std::vector<float> &w;  // directional weights
    const std::array<float, 3> &extf;  // external force

    float omega;  // parameters for CollisionSRT dynamics
    float tau;

    std::array<float, 3> _get_updated_u(
        const size_t zl, const size_t yl, const size_t xl,
        const float rholocal, const float *ulocal) const;
    void _get_cu(const std::array<float, 3> &u, float *cu) const;
    void _collision_ref(SimData &simdata) const;
    void _collision_tbb(SimData &simdata) const;
    void _collision_kernel_avx2(
        const size_t zl, const size_t yl, const size_t xl,
        SimData &simdata,
        float *scratch) const;  // scratch space to compute dot(ck,u)
    void _collision_kernel(
        const size_t zl, const size_t yl, const size_t xl,
        SimData &simdata,
        float *scratch) const;  // scratch space to compute dot(ck,u)

 public:
    CollisionSRT() = delete;
    CollisionSRT(const LBModel *lbmodel, const Domain *domain);
    CollisionSRT(CollisionSRT&) = delete;
    CollisionSRT& operator=(CollisionSRT&) = delete;
    ~CollisionSRT();
    void operator()(SimData &simdata, bool reference = false) const;
};

#endif
