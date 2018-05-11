#ifndef BGK_HPP
#define BGK_HPP

#include "LBDynamics.hpp"
#include "LBModel.hpp"
#include "Domain.hpp"

class BGK: public LBDynamics{
    
private:
    
    // It's convenient to store the input const references
    const LBModel *_lbmodel;
    const Domain *_domain;

    float _omega;  // parameters for BGK dynamics
    float _tau;

    void _collide_ref(Lattice &lattice) const;
    void _collide(Lattice &lattice) const;
    void _collide_kernel(
        const size_t zyx,
        const size_t kdim,
        const std::vector<int32_t> &c, // lattice velocities
        const std::vector<float> &w, // directional weights
        const std::array<float, 3> &extForce,
        float * __restrict__ n,
        const float * __restrict__ rho,
        const float * __restrict__ u,
        float *cu // scratch space to compute dot(ck,u)
    ) const;
    void _stream_ref(Lattice &lattice) const;
    void _stream(Lattice &lattice) const;

public:
    
    BGK() = delete;
    BGK(const LBModel *lbmodel, const Domain *domain);
    ~BGK();
    void collideAndStream(Lattice &lattice) const;
    void calcMoments(Lattice &lattice) const;

};

#endif
