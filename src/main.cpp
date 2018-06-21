#include <iostream>
#include <memory>
#include <iomanip>

#include "tbb/tbb.h"

#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "Flow.hpp"
#include "SimData.hpp"
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
    
    // Simulation data
    const size_t kdim = lbmodel->get_num_directions();
    auto simdata = SimData(domain->get_dimensions(), kdim);

    const auto calc_moments = std::make_unique<CalcMoments>();

    // Initialize
    domain->initialize(simdata);
    boundary->reset(simdata);
    simdata.write_state("InitState.h5");

    // Time loop
    const tbb::tick_count start = tbb::tick_count::now();
    for (auto istep=0; istep<num_timesteps; ++istep){
        collide->operator()(simdata);
        stream->operator()(simdata);
        boundary->apply_noslip(simdata);
        boundary->apply_periodicity(simdata);
        calc_moments->operator()(lbmodel, simdata);
        boundary->reset(simdata);
    }
    const auto elapsed = (tbb::tick_count::now()-start).seconds();
    std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
    simdata.write_state("FinalState.h5");

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
