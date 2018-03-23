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

    auto flow = std::make_unique<Flow>("Couette");
    auto domain = flow->getFlowDomain();
    auto boundary = flow->getFlowBoundary();
    auto lbmodel = LBModel("D3Q27");
    auto lattice = Lattice(lbmodel, domain);
    boundary->apply(lattice);
    auto lbdynamics = std::make_unique<BGK>(lbmodel, domain);
    
    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for(auto i=0; i<100; ++i){
        lbdynamics->collideAndStream(lattice);
        boundary->apply(lattice);
    }
    std::cout<<"Time: "<<(tbb::tick_count::now()-start).seconds()<<"s\n";
    lattice.writeState();
    
    // Finalize
    
    return 0;

}
