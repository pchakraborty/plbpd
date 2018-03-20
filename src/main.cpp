#include <iostream>
#include <memory>
#include "LBModel.hpp"
#include "Domain.hpp"
#include "Lattice.hpp"
#include "LBDynamics.hpp"
#include "BGK.hpp"
#include "tbb/tbb.h"

int main(){

    auto nthreads = 4; // tbb::task_scheduler_init::default_num_threads();
    tbb::task_scheduler_init init(nthreads);
    std::cout<<"plbpd: using "<<nthreads<<" threads..."<<std::endl;
    
    auto lbmodel = LBModel("D3Q27");
    auto domain = Domain("domain.yml");
    auto lattice = Lattice(lbmodel, domain, true); // bootstrap = true
    auto lbdynamics = std::make_unique<BGK>(lbmodel, domain, lattice);
    
    // Time loop
    tbb::tick_count start = tbb::tick_count::now();
    for(auto i=0; i<100; ++i)
        lbdynamics->collideAndStream();
    std::cout<<"Time: "<<(tbb::tick_count::now()-start).seconds()<<"s\n";
    // Finalize
    
    return 0;

}
