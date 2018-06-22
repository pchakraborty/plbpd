#include <iostream>
#include <memory>
#include <iomanip>

#include "tbb/tbb.h"

#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "Flow.hpp"
#include "SimData.hpp"
#include "LBDynamics.hpp"
#include "CalcMoments.hpp"

void initialize(const Flow* flow, SimData& simdata){
    const auto domain = flow->get_domain();
    const auto boundary = flow->get_boundary();
    domain->initialize(simdata);
    boundary->reset(simdata);
}

void run_timeloop(
    const Flow* flow,
    const LBDynamics* lbdynamics,
    const CalcMoments* calc_moments,
    SimData& simdata){
    // Timeloop
    const auto boundary = flow->get_boundary();
    const auto lbmodel = flow->get_lbmodel();
    auto num_timesteps = flow->get_num_timesteps();
    for (auto istep=0; istep<num_timesteps; ++istep){
        lbdynamics->collide(simdata);
        lbdynamics->stream(simdata);
        boundary->apply_noslip(simdata);
        boundary->apply_periodicity(simdata);
        calc_moments->operator()(lbmodel, simdata);
        boundary->reset(simdata);
    }
}

void finalize(const Flow* flow, const LBDynamics* lbdynamics, const CalcMoments* calc_moments){
    auto boundary = flow->get_boundary();
    // Print times
    std::cout<<"LBDynamics: "<<lbdynamics->get_total_time()<<"s\n";
    std::cout<<"-collide: "<<lbdynamics->get_time_collide()<<"s\n";
    std::cout<<"-stream: "<<lbdynamics->get_time_stream()<<"s\n";
    std::cout<<"CalcMoments: "<<calc_moments->get_total_time()<<"s\n";
    std::cout<<"Boundary: "<<boundary->get_total_time()<<"s\n";
    std::cout<<"-noslip: "<<boundary->get_time_noslip()<<"s\n";
    std::cout<<"-periodicity: "<<boundary->get_time_periodicity()<<"s\n";
    std::cout<<"-reset: "<<boundary->get_time_reset()<<"s\n";
}    

int main(){

    const auto flow = std::make_unique<Flow>("Poiseuille2D");
    const auto domain = flow->get_domain();
    const auto lbmodel = flow->get_lbmodel();

    const auto lbdynamics = std::make_unique<LBDynamics>(lbmodel, domain);
    auto simdata = SimData(domain->get_dimensions(), lbmodel->get_num_directions());
    const auto calc_moments = std::make_unique<CalcMoments>();

    initialize(flow.get(), simdata);
    {
        simdata.write_state("InitState.h5");
        const tbb::tick_count start = tbb::tick_count::now();
        run_timeloop(flow.get(), lbdynamics.get(), calc_moments.get(), simdata);
        const auto elapsed = (tbb::tick_count::now()-start).seconds();
        std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
        simdata.write_state("FinalState.h5");
    }
    finalize(flow.get(), lbdynamics.get(), calc_moments.get());

    return 0;

}
