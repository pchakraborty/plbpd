#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "LBModel.hpp"
// #include "ArrayND.hpp"
#include "Field.hpp"

class Lattice final{

private:

    const std::tuple<size_t, size_t, size_t> _domain_dimensions;
    const size_t _kdim; // LBModel::num_directions

    void _bootstrap();
    void _read_state();
    
public:

    // array3f *rho; // density at each node
    // array4f *u; // velocity, u[3] at each node
    // array4f *n; // particle distribution, n[num_directions] at each node
    // array4f *ntmp; // for streaming

    Field::ScalarField<float, 1> *rho; // 1 -> number of buffer layers
    Field::VectorField<float, 1> *u;
    Field::VectorField<float, 1> *n;
    Field::VectorField<float, 1> *ntmp;
    
    Lattice() = delete;
    Lattice(
        const std::tuple<size_t, size_t, size_t> domain_dimensions, 
        const size_t num_directions
    );
    ~Lattice();
    void write_state(std::string dump_file);

};

#endif
