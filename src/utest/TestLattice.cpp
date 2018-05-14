#include "../Lattice.hpp"
#include <tuple>
#include <memory>

int main(){

    const size_t zdim=5, ydim=4, xdim=3, kdim=19;
    auto zyxdims = std::make_tuple(zdim, ydim, xdim);
    auto lattice = std::make_unique<Lattice>(zyxdims, kdim);
    lattice->writeState("BootstrappedState.h5");
    
    return 0;
}
