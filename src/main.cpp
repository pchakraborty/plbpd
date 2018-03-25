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
    auto flow = std::make_unique<Flow>("Couette");
    auto domain = flow->getFlowDomain();
    auto boundary = flow->getFlowBoundary();

    // Lattice Boltzmann model
    auto lbmodel = LBModel("D3Q27");

    // Lattice Boltzmann dynamics
    auto lbdynamics = std::make_unique<BGK>(lbmodel, domain);
    
    // Initialize data
    auto lattice = Lattice(lbmodel, domain);
    lbdynamics->initialize(lattice);
    boundary->apply(domain, lattice);

    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for(auto i=0; i<100; ++i){
        lbdynamics->collideAndStream(lattice);
        boundary->apply(domain, lattice);
    }
    std::cout<<"Time: "<<(tbb::tick_count::now()-start).seconds()<<"s\n";
    lattice.writeState();
    
    // Finalize
    
    return 0;

}
