#ifndef SIMDATA_HPP
#define SIMDATA_HPP

#include "Field.hpp"

class SimData final{

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
    
    SimData() = delete;
    SimData(
        const std::tuple<size_t, size_t, size_t> domain_dimensions, 
        const size_t num_directions
    );
    ~SimData();
    void write_state(std::string dump_file);

};

#endif
