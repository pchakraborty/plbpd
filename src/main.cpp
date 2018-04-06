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

    // Lattice Boltzmann model
    auto lbmodel = std::make_unique<LBModel>("D3Q27");
    
    // Lattice Boltzmann dynamics
    auto lbdynamics = std::make_unique<BGK>(lbmodel.get(), domain);
    
    // Initialize problem
    auto lattice = Lattice(lbmodel.get(), domain);
    lbdynamics->initialize(lattice);
    boundary->apply(lattice);

    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for (auto i=0; i<100; ++i){
        lbdynamics->collideAndStream(lattice);
        // lbdynamics->calcMoments(lattice);
        boundary->apply(lattice);
    }
    std::cout<<"Time: "<<(tbb::tick_count::now()-start).seconds()<<"s"<<std::endl;
    lattice.writeState();
    
    // Finalize
    return 0;

}
