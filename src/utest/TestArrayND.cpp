#include <iostream>
#include "../ArrayND.hpp"
#undef NDEBUG
#include <cassert>

int main(){

    const auto zdim = 5;
    const auto ydim = 4;
    const auto xdim = 3;
    const auto kdim = 2;

    // 3D array
    ArrayND::Array3D<float> A3(zdim, ydim, xdim, 0.0f);
    for (auto zl=0; zl<zdim; ++zl)
        for (auto yl=0; yl<ydim; ++yl)
            for (auto xl=0; xl<xdim; ++xl)
                A3.at(zl,yl,xl) = xl+(yl+zl*ydim)*xdim;
    std::vector<float> A3v(zdim*ydim*xdim, 0.0f);
    for (auto i=0; i!=A3v.size(); ++i)
        A3v[i] = i;
    float error3 = 0.0f;
    for (auto zl=0; zl<zdim; ++zl)
        for (auto yl=0; yl<ydim; ++yl)
            for (auto xl=0; xl<xdim; ++xl)
                error3 += A3.at(zl,yl,xl) - A3v[xl+(yl+zl*ydim)*xdim];
    assert(error3==0.0f);

    // 4D array
    ArrayND::Array4D<float> A4(zdim, ydim, xdim,kdim, 0.0f);
    for (auto zl=0; zl<zdim; ++zl)
        for (auto yl=0; yl<ydim; ++yl)
            for (auto xl=0; xl<xdim; ++xl){
                auto ndx3 = xl+(yl+zl*ydim)*xdim;
                for (auto kl=0; kl<kdim; ++kl)
                    A4.at(zl,yl,xl,kl) = ndx3*kdim + kl;
            }
    std::vector<float> A4v(zdim*ydim*xdim*kdim, 0.0f);
    for (auto i=0; i!=A4v.size(); ++i)
        A4v[i] = i;
    float error4 = 0.0f;
    for (auto zl=0; zl<zdim; ++zl)
        for (auto yl=0; yl<ydim; ++yl)
            for (auto xl=0; xl<xdim; ++xl){
                auto ndx3 = xl+(yl+zl*ydim)*xdim;
                for (auto kl=0; kl<kdim; ++kl)
                    error4 += A4.at(zl,yl,xl,kl) - A4v[ndx3*kdim+kl];
            }
    assert(error4==0.0f);
    
    return 0;
}
