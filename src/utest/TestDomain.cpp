#include "../Domain.hpp"
#include "../LBModel.hpp"
#include <memory>
#include <cassert>


int main(){
    const auto lbmodel = std::make_unique<LBModel>("D3Q19");

    const size_t xdim = 5, ydim = 4, zdim = 3;
    const float nu = 0.5, rho = 1.0;
    const std::array<float, 3> vel = {-1148.14, 762.878, 3.2};
    const std::array<float, 3> force = {-0.145, 0.0, 4.2};
    const auto domain = std::make_unique<Domain>
        (lbmodel.get(), xdim, ydim, zdim, nu, rho, vel, force);
    assert(std::tie(xdim, ydim, zdim)==domain->getDimensions());
    assert(nu=domain->getFluidViscosity());
    assert(rho=domain->getFluidDensity());
    assert(vel==domain->getInitFlowVelocity());
    assert(force==domain->getExternalForce());
    
    return 0;
}
    
