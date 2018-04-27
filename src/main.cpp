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

    // Flow details: domain, boundary, model
    auto flow = std::make_unique<Flow>("Couette3D");
    auto domain = flow->getFlowDomain();
    auto boundary = flow->getFlowBoundary();
    auto lbmodel = flow->getLBModel();
    auto numTimeSteps = flow->getNumTimeSteps();
    
    // Lattice Boltzmann dynamics
    auto lbdynamics = std::make_unique<BGK>(lbmodel, domain);
    
    // Lattice
    size_t kdim = lbmodel->getNumberOfDirections();
    auto lattice = Lattice(domain->getDimensions(), kdim);

    // Initialize
    domain->initialize(lattice);
    boundary->reset(lattice);
    lbdynamics->initialize(lattice);
    lbdynamics->calcMoments(lattice);
    lattice.writeState("InitState.h5");
    
    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for (auto i=0; i<numTimeSteps; ++i){
        // std::cout<<"step: "<<i<<std::endl;
        lbdynamics->collideAndStream(lattice);
        boundary->applyNoslip(lattice);
        boundary->applyPeriodicity(lattice);
        lbdynamics->calcMoments(lattice);
        boundary->reset(lattice);
    }
    auto elapsed = (tbb::tick_count::now()-start).seconds();
    std::cout<<"Time taken by timeloop: "<<elapsed<<"s"<<std::endl;
    lattice.writeState("FinalState.h5");
    
    // Finalize
    return 0;

}
