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
    Lattice &_lattice;
    
    float omega;  // parameters for BGK dynamics
    float tau;
    
    //__attribute__((always_inline))
    std::array<float, 3> _getEqlbVelocity(
      size_t zl, size_t yl, size_t xl,
      size_t zdim, size_t ydim, size_t xdim
    );
    void _printInfoForDebugging();
    void _collideAndStreamOnPlane(size_t zl);
    void _serialCollideAndStream();
    void _parallelCollideAndStream();
    
public:
    
    BGK() = delete;
    BGK(const LBModel &lbmodel, const Domain &domain, Lattice &lattice);
    ~BGK();
    void setup();
    void collideAndStream();
    void calcMoments();
    float getAvgFluidDensity();
    void getEqlbDist(const float rholoc, const float &uloc, float &nloc);
    
};

#endif
