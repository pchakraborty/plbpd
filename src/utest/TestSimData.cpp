#include <tuple>
#include <memory>

#include "../SimData.hpp"

#undef NDEBUG

int main() {
    const size_t zdim = 5, ydim = 4, xdim = 3, kdim = 19;
    auto zyxdims = std::make_tuple(zdim, ydim, xdim);
    auto simdata = std::make_unique<SimData>(zyxdims, kdim);
    simdata->write_state("BootstrappedState.h5");

    return 0;
}
