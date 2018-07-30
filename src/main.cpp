#include <iostream>
#include <memory>
#include <iomanip>

#include "tbb/tbb.h"

#include "LBModel.hpp"
#include "Domain.hpp"
#include "Patch.hpp"
#include "Boundary.hpp"
#include "Flow.hpp"
#include "SimData.hpp"
#include "LBDynamics.hpp"
#include "CalcMoments.hpp"

using FlowDynamicsMoments = std::tuple<const Flow*, const LBDynamics*, const CalcMoments*>;

void initialize(const Flow* flow, SimData& simdata){
    const auto domain = flow->get_domain();
    const auto boundary = flow->get_boundary();
    domain->initialize(simdata);
    boundary->reset(simdata);
}

void run_timeloop(FlowDynamicsMoments flow_dyn_mom, SimData& simdata){
    auto flow = std::get<0>(flow_dyn_mom);
    auto lbdynamics = std::get<1>(flow_dyn_mom);
    auto calc_moments = std::get<2>(flow_dyn_mom);
    const auto boundary = flow->get_boundary();
    const auto lbmodel = flow->get_lbmodel();
    auto num_timesteps = flow->get_num_timesteps();
    // Timeloop
    for (auto istep=0; istep<num_timesteps; ++istep){
        lbdynamics->collide(simdata);
        lbdynamics->stream(simdata);
        boundary->apply_noslip(simdata);
        boundary->apply_periodicity(simdata);
        calc_moments->operator()(lbmodel, simdata);
        boundary->reset(simdata);
    }
}

void finalize(FlowDynamicsMoments flow_dyn_mom){
    const auto flow = std::get<0>(flow_dyn_mom);
    const auto lbdynamics = std::get<1>(flow_dyn_mom);
    const auto calc_moments = std::get<2>(flow_dyn_mom);
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
    auto flow_dyn_mom = std::make_tuple(flow.get(), lbdynamics.get(), calc_moments.get());

    simdata.write_state("InitState.h5");

    const tbb::tick_count start = tbb::tick_count::now();
    run_timeloop(flow_dyn_mom, simdata);
    const auto elapsed = (tbb::tick_count::now()-start).seconds();

    simdata.write_state("FinalState.h5");

    std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
    finalize(flow_dyn_mom);

    return 0;

}
