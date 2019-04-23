#include "../Block3d.hpp"

#undef NDEBUG

int main() {
    const size_t xsize = 5, ysize = 4, zsize = 3;
    const auto domain = std::make_unique<Domain>
        (lbmodel.get(), xdim, ydim, zdim, nu, rho, vel, force);
    assert(std::tie(xdim, ydim, zdim) == domain->get_dimensions());
    return 0;
}
