#include <iostream>
#include <memory>
#include <iomanip>

#include "tbb/tbb.h"

#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "Flow.hpp"
#include "Lattice.hpp"
#include "CollisionSRT.hpp"
#include "Streaming.hpp"
#include "CalcMoments.hpp"

int main(){

    // Flow definition: domain, boundary, model
    const auto flow = std::make_unique<Flow>("Poiseuille2D");
    const auto domain = flow->get_flow_domain();
    const auto boundary = flow->get_flow_boundary();
    const auto lbmodel = flow->get_lbmodel();
    const auto num_timesteps = flow->get_num_timesteps();

    // Lattice Boltzmann dynamics
    const auto collide = std::make_unique<CollisionSRT>(lbmodel, domain);
    const auto stream = std::make_unique<Streaming>(lbmodel, "push");
    
    // Lattice
    const size_t kdim = lbmodel->get_num_directions();
    auto lattice = Lattice(domain->get_dimensions(), kdim);

    const auto calc_moments = std::make_unique<CalcMoments>();

    // Initialize
    domain->initialize(lattice);
    boundary->reset(lattice);
    lattice.write_state("InitState.h5");

    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for (auto i=0; i<num_timesteps; ++i){
        collide->operator()(lattice);
        stream->operator()(lattice);
        boundary->apply_noslip(lattice);
        boundary->apply_periodicity(lattice);
        calc_moments->operator()(lbmodel, lattice);
        boundary->reset(lattice);
    }
    auto elapsed = (tbb::tick_count::now()-start).seconds();
    std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
    lattice.write_state("FinalState.h5");

    // Print times
    std::cout<<"Collide: "<<collide->get_total_time()<<"s\n";
    std::cout<<"Stream: "<<stream->get_total_time()<<"s\n";
    std::cout<<"CalcMoments: "<<calc_moments->get_total_time()<<"s\n";
    std::cout<<"Boundary: "<<boundary->get_total_time()<<"s\n";
    std::cout<<"-noslip: "<<boundary->get_time_noslip()<<"s\n";
    std::cout<<"-periodicity: "<<boundary->get_time_periodicity()<<"s\n";
    std::cout<<"-reset: "<<boundary->get_time_reset()<<"s\n";

    // Finalize
    return 0;

}
