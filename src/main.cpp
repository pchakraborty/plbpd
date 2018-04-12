#include <iostream>
#include <memory>
#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "Flow.hpp"
#include "Lattice.hpp"
#include "LBDynamics.hpp"
#include "BGK.hpp"
#include "tbb/tbb.h"

int main(){

    // Domain and its boundary
    // Domain is a collection of connected patches
    auto flow = std::make_unique<Flow>("Couette");
    auto domain = flow->getFlowDomain();
    auto boundary = flow->getFlowBoundary();
    auto lbmodel = flow->getLBModel();
    
    // Lattice Boltzmann dynamics
    auto lbdynamics = std::make_unique<BGK>(lbmodel, domain);
    
    // Initialize problem
    auto lattice = Lattice(lbmodel, domain);
    lbdynamics->initialize(lattice);
    boundary->apply(lattice);
    lattice.writeState("InitState.h5");
    
    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for (auto i=0; i<100; ++i){
        lbdynamics->collideAndStream(lattice);
        lbdynamics->calcMoments(lattice);
        boundary->apply(lattice);
    }
    auto elapsed = (tbb::tick_count::now()-start).seconds();
    std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
    lattice.writeState("FinalState.h5");
    
    // Finalize
    return 0;

}
