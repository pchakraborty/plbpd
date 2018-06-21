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
    auto flow = std::make_unique<Flow>("Poiseuille2D");
    auto domain = flow->get_flow_domain();
    auto boundary = flow->get_flow_boundary();
    auto lbmodel = flow->get_lbmodel();
    auto num_timesteps = flow->get_num_timesteps();

    // Lattice Boltzmann dynamics
    auto collide = std::make_unique<CollisionSRT>(lbmodel, domain);
    auto stream = std::make_unique<Streaming>(lbmodel, "push");
    
    // Lattice
    size_t kdim = lbmodel->get_num_directions();
    auto lattice = Lattice(domain->get_dimensions(), kdim);

    auto calc_moments = CalcMoments();

    // Initialize
    domain->initialize(lattice);
    boundary->reset(lattice);
    lattice.write_state("InitState.h5");

    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for (auto i=0; i<num_timesteps; ++i){
        // std::cout<<"step: "<<i<<std::endl;
        (*collide)(lattice);
        (*stream)(lattice);
        boundary->apply_noslip(lattice);
        boundary->apply_periodicity(lattice);
        calc_moments(lbmodel, lattice);
        boundary->reset(lattice);
    }
    auto elapsed = (tbb::tick_count::now()-start).seconds();
    std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
    lattice.write_state("FinalState.h5");

    // Print times
    std::cout<<"Collide: "<<collide->get_total_time()<<"s\n";
    std::cout<<"Stream: "<<stream->get_total_time()<<"s\n";
    std::cout<<"CalcMoments: "<<calc_moments.get_total_time()<<"s\n";
    std::cout<<"Boundary: "<<boundary->get_total_time()<<"s\n";
    std::cout<<"-noslip: "<<boundary->get_time_noslip()<<"s\n";
    std::cout<<"-periodicity: "<<boundary->get_time_periodicity()<<"s\n";
    std::cout<<"-reset: "<<boundary->get_time_reset()<<"s\n";

    // Finalize
    return 0;

}
