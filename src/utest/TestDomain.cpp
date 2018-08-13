#include <memory>
#include <cassert>

#include "../Domain.hpp"
#include "../LBModel.hpp"

#undef NDEBUG

int main() {
    const auto lbmodel = std::make_unique<LBModel>("D3Q19");

    const size_t xdim = 5, ydim = 4, zdim = 3;
    const float nu = 0.5, rho = 1.0;
    const std::array<float, 3> vel = {-1148.14, 762.878, 3.2};
    const std::array<float, 3> force = {-0.145, 0.0, 4.2};
    const auto domain = std::make_unique<Domain>
        (lbmodel.get(), xdim, ydim, zdim, nu, rho, vel, force);
    assert(std::tie(xdim, ydim, zdim) == domain->get_dimensions());
    assert(nu == domain->get_fluid_viscosity());
    assert(rho == domain->get_fluid_density());
    assert(vel == domain->get_init_flow_velocity());
    assert(force == domain->get_external_force());
    return 0;
}
