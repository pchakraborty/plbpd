#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "LBModel.hpp"
#include "Domain.hpp"

#include "ArrayND.hpp"
using array3f = ArrayND::Array3D<float>;
using array4f = ArrayND::Array4D_simd<float>;

class Lattice final{

private:

    // It's convenient to store input const references
    const LBModel *_lbmodel;
    const Domain *_domain;
    
    void _bootstrap();
    void _readState();
    
public:

    array3f *rho; // density at each node
    array4f *u; // velocity, u[3] at each node
    array4f *n; // particle distribution, n[numberOfDirections] at each node
    array4f *ntmp; // for streaming
        
    Lattice() = delete;
    Lattice(const LBModel *lbmodel, const Domain *domain);
    ~Lattice();
    void writeState(std::string dumpFile);

};

#endif
