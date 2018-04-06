#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "LBModel.hpp"
#include "Domain.hpp"
#include <memory>

class Lattice final{

private:

    // It's convenient to store input const references
    const LBModel *_lbmodel;
    const Domain *_domain;
    
    void _bootstrap();
    void _readState();
    
public:

    // Data storage
    float *rho; // density at each node
    float *u; // velocity, u[3], at each node
    float *n; // particle distribution, n[numVelocityVectors] at each node
    float *ntmp; // for streaming
        
    Lattice() = delete;
    Lattice(const LBModel *lbmodel, const Domain *domain);
    ~Lattice();
    void writeState();

};

#endif
