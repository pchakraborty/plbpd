#ifndef COLLISION_SRT_HPP
#define COLLISION_SRT_HPP

#include "Collision.hpp"
#include "LBModel.hpp"
#include "Domain.hpp"

class CollisionSRT final: public Collision{
    
private:
    
    // It's convenient to store the input const references
    const LBModel *_lbmodel;
    const Domain *_domain;

    bool _reference;
    
    float _omega;  // parameters for CollisionSRT dynamics
    float _tau;

    void _collision_ref(SimData &simdata) const;
    void _collision_tbb(SimData &simdata) const;
    void _collision_kernel_avx2(
        const size_t zyx,
        const size_t kdim,
        const std::vector<int32_t> &c, // directional velocities
        const std::vector<float> &w, // directional weights
        const std::array<float, 3> &ext_force,
        float * __restrict__ n,
        const float * __restrict__ rho,
        const float * __restrict__ u,
        float *cu // scratch space to compute dot(ck,u)
    ) const;
    void _collision_kernel(
        const size_t zyx,
        const size_t kdim,
        const std::vector<int32_t> &c, // directional velocities
        const std::vector<float> &w, // directional weights
        const std::array<float, 3> &extForce,
        float * __restrict__ n,
        const float * __restrict__ rho,
        const float * __restrict__ u,
        float *cu // scratch space to compute dot(ck,u)
    ) const;

public:
    
    CollisionSRT() = delete;
    CollisionSRT(const LBModel *lbmodel, const Domain *domain);
    CollisionSRT(const LBModel *lbmodel, const Domain *domain, bool reference);
    CollisionSRT(CollisionSRT&) = delete;
    CollisionSRT& operator=(CollisionSRT&) = delete;
    ~CollisionSRT();
    void operator()(SimData &simdata) const;

};

#endif
