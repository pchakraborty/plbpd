#include <iostream>
#include <memory>
#include <iomanip>

#include "tbb/tbb.h"

#include "Flow.hpp"
#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "SimData.hpp"
#include "Plbpd.hpp"

int main(){
    // Model, Domain and Boundary
    const auto flow = std::make_unique<Flow>("Poiseuille2D");
    const auto lbmodel = flow->get_lbmodel();
    const auto domain = flow->get_domain();
    const auto boundary = flow->get_boundary();

    // Define problem and its data
    auto my_problem = std::make_unique<Plbpd>(lbmodel, domain, boundary);
    auto simdata = SimData(domain->get_dimensions(), lbmodel->get_num_directions());
    
    // Solve problem
    my_problem->initialize(simdata);
    simdata.write_state("init_state.h5");
    my_problem->march_in_time(flow->get_num_timesteps(), simdata);
    simdata.write_state("final_state.h5");
    my_problem->print_times();
}
