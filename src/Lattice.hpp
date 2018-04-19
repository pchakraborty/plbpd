#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "LBModel.hpp"
#include "ArrayND.hpp"
using array3f = ArrayND::Array3D<float>;
using array4f = ArrayND::Array4D<float>;

class Lattice final{

private:

    const std::tuple<size_t, size_t, size_t> _domainDimensions;
    const size_t _kdim; // LBModel::numberOfDirections

    void _bootstrap();
    void _readState();
    
public:

    array3f *rho; // density at each node
    array4f *u; // velocity, u[3] at each node
    array4f *n; // particle distribution, n[numberOfDirections] at each node
    array4f *ntmp; // for streaming
        
    Lattice() = delete;
    Lattice(
        const std::tuple<size_t, size_t, size_t> domainDimensions, 
        const size_t numberOfDirections
    );
    ~Lattice();
    void writeState(std::string dumpFile);

};

#endif
