#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "LBModel.hpp"
#include "Domain.hpp"

class Lattice final{

private:

    // It's convenient to store input const references
    const LBModel *_lbmodel;
    const Domain *_domain;
    
    void _bootstrap();
    void _readState();
    
public:

    // Data storage
    std::vector<unsigned short> nodetype; // one of {BUFFER, INTFLD, etc.}
    std::vector<float> rho; // density at each node
    std::vector<float> u; // velocity, u[3], at each node
    std::vector<float> n; // particle distribution, n[numVelocityVectors] at each node
    std::vector<float> ntmp; // for streaming
    
    Lattice() = delete;
    Lattice(const LBModel *lbmodel, const Domain *domain);
    ~Lattice();
    void writeState();

};

#endif
