#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <memory>
#include "Domain.hpp"
#include "LBModel.hpp"

class Lattice final{

private:

    size_t _numVelocityVectors; // local copy of LBModel::_numVelocityVectors
    size_t _xdim, _ydim, _zdim; // local copies of Domain dimensions
    
    void _bootstrapRestarts();
    void _readRestarts(std::string restartFile);
    
public:

    // Data storage
    std::vector<unsigned short> nodetype; // one of {BUFFER, INTFLD, etc.}
    std::vector<double> rho; // density at each node
    std::vector<double> u; // velocity, u[3], at each node
    std::vector<double> n; // particle distribution, n[numVelocityVectors] at each node
    std::vector<double> ntmp; // for streaming
    
    Lattice() = delete;
    Lattice(const LBModel &lbmodel, const Domain &domain, bool bootstrap);
    ~Lattice();
    void dumpState();
    
    
};

#endif
