#include <iostream>
#include <memory>
#include <iomanip>

#include "tbb/tbb.h"

#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "Flow.hpp"
#include "Lattice.hpp"
#include "LBDynamics.hpp"
#include "BGK.hpp"
#include "CalcMoments.hpp"

int main(){

    // Flow definition: domain, boundary, model
    auto flow = std::make_unique<Flow>("Poiseuille2D");
    auto domain = flow->get_flow_domain();
    auto boundary = flow->get_flow_boundary();
    auto lbmodel = flow->get_lbmodel();
    auto num_timesteps = flow->get_num_timesteps();

    // Lattice Boltzmann dynamics
    auto lbdynamics = std::make_unique<BGK>(lbmodel, domain);

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
        lbdynamics->collide(lattice);
        lbdynamics->stream(lattice);
        boundary->apply_noslip(lattice);
        boundary->apply_periodicity(lattice);
        calc_moments(lbmodel, lattice);
        boundary->reset(lattice);
    }
    auto elapsed = (tbb::tick_count::now()-start).seconds();
    std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
    lattice.write_state("FinalState.h5");

    // Print times
    std::cout<<"LBDynamics: "<<lbdynamics->get_total_time()<<"s\n";
    std::cout<<"-collide: "<<lbdynamics->get_time_collide()<<"s\n";
    std::cout<<"-stream: "<<lbdynamics->get_time_stream()<<"s\n";
    std::cout<<"CalcMoments: "<<calc_moments.get_total_time()<<"s\n";
    std::cout<<"Boundary: "<<boundary->get_total_time()<<"s\n";
    std::cout<<"-noslip: "<<boundary->get_time_noslip()<<"s\n";
    std::cout<<"-periodicity: "<<boundary->get_time_periodicity()<<"s\n";
    std::cout<<"-reset: "<<boundary->get_time_reset()<<"s\n";

    // Finalize
    return 0;

}
