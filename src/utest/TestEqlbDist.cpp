#include <iostream>
#include <memory>
#undef NDEBUG
#include <cassert>
#include "../LBModel.hpp"
#include "../EqlbDist.hpp"

int main(){
    auto lbmodel = std::make_unique<LBModel>("D2Q9");
    float rholocal = 1.0f;
    std::array<float, 3> ulocal = {0.5, 0.5, -0.5};
    std::vector<float> nlocal(lbmodel->get_num_directions(), 0.0f);

    // Get eqlb dist
    auto eqlbdist = std::make_unique<EqlbDist>();
    (*eqlbdist)(lbmodel.get(), rholocal, ulocal, nlocal);

    // Expected value of eqlb dist
    std::vector<float> nexpected = {
        -0.05555555597, 0.277777791, -0.05555555597, -0.05555555597, 0.277777791,
        -0.003472222248, -0.003472222248, 0.2048611194, 0.03819444403
    };

    // Error
    auto err = 0.0f;
    for (auto k=0; k<9; ++k){
        auto tmp = nexpected[k]-nlocal[k];
        err += tmp*tmp;
    }
    assert(err==0.0f);
}
