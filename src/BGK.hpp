#ifndef BGK_HPP
#define BGK_HPP

#include "LBDynamics.hpp"
#include "LBModel.hpp"
#include "Lattice.hpp"

class BGK: public LBDynamics{
    
private:
    
    // It's convenient to store the input const references
    const LBModel &_lbmodel;
    const Domain &_domain;
    
    float omega;  // parameters for BGK dynamics
    float tau;
    
    //__attribute__((always_inline))
    std::array<float, 3> _getEqlbVelocity(
      size_t zl, size_t yl, size_t xl,
      size_t zdim, size_t ydim, size_t xdim,
      Lattice &lattice
    );
    void _printInfoForDebugging();
    void _collideAndStreamOnPlane(size_t zl, Lattice &lattice);
    void _serialCollideAndStream(Lattice &lattice);
    void _parallelCollideAndStream(Lattice &lattice);
    
public:
    
    BGK() = delete;
    BGK(const LBModel &lbmodel, const Domain &domain);
    ~BGK();
    void collideAndStream(Lattice &lattice);
    void calcMoments();
    float getAvgFluidDensity();
    void getEqlbDist(const float rholoc, const float &uloc, float &nloc);
    
};

#endif
