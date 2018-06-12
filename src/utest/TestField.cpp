#include <iostream>
#include "../Field.hpp"
#undef NDEBUG
#include <cassert>
#include <vector>
#include <memory>

int main(){
    
    {
        // Scalar field with 0 buffer layers

        const auto zsize = 4;
        const auto ysize = 5;
        const auto xsize = 7;
        
        // Scalar field
        Field::Field<float,0> sf(zsize, ysize, xsize, 0.0f);
        size_t zdim, ydim, xdim;
        std::tie(zdim, ydim, xdim) = sf.get_dimensions();
        assert(std::tie(zdim, ydim, xdim) == std::tie(zsize, ysize, xsize));
        assert(sf.get_vector_length() == 1);
        // Set data using iterator
        auto ctr = 0;
        for (auto it=sf.begin(); it!=sf.end(); ++it)
            *it = ctr++;
        // Get data using at()
        Field::FieldExtents e = sf.get_extents();
        assert(e.zbegin == 0); assert(e.zend == zsize);
        assert(e.ybegin == 0); assert(e.yend == ysize);
        assert(e.xbegin == 0); assert(e.xend == xsize);
        for (auto zl=e.zbegin; zl<e.zend; ++zl)
            for (auto yl=e.ybegin; yl<e.yend; ++yl)
                for (auto xl=e.xbegin; xl<e.xend; ++xl){
                    auto ndx3d = xl+(yl+zl*ydim)*xdim;
                    assert(sf.at(zl,yl,xl) == ndx3d);
                }
    }
    
    {
        // Vector field with 0 buffer layers

        const auto zsize = 4;
        const auto ysize = 5;
        const auto xsize = 7;
        const auto vsize = 2;

        Field::Field<float,0> vf(zsize, ysize, xsize, vsize, 0.0f);
        size_t zdim, ydim, xdim;
        std::tie(zdim, ydim, xdim) = vf.get_dimensions();
        assert(std::tie(zdim, ydim, xdim) == std::tie(zsize, ysize, xsize));
        assert(vf.get_vector_length() == vsize);
        // Set data using iterator
        auto ctr = 0;
        for (auto it=vf.begin(); it!=vf.end(); ++it)
            *it = ctr++;

        // Get data using at()
        Field::FieldExtents e = vf.get_extents();
        assert(e.zbegin == 0); assert(e.zend == zsize);
        assert(e.ybegin == 0); assert(e.yend == ysize);
        assert(e.xbegin == 0); assert(e.xend == xsize);
        auto vlen = vf.get_vector_length();
        for (auto zl=e.zbegin; zl<e.zend; ++zl)
            for (auto yl=e.ybegin; yl<e.yend; ++yl)
                for (auto xl=e.xbegin; xl<e.xend; ++xl){
                    auto ndx3d = xl+(yl+zl*ydim)*xdim;
                    for (auto kl=0; kl<vlen; ++kl)
                        assert(vf.at(zl,yl,xl,kl) == kl+ndx3d*vlen);
                }
    }

    {
        const auto zsize=2;
        const auto ysize=2;
        const auto xsize=3;
        
        // Scalar field with one buffer layer
        auto sf = std::make_unique<Field::Field<float,1> >(zsize, ysize, xsize, 0.0f);
        size_t zdim, ydim, xdim;
        std::tie(zdim, ydim, xdim) = sf->get_dimensions();
        assert(zdim == zsize+2); assert(ydim == ysize+2); assert(xdim == xsize+2);
        Field::FieldExtents e = sf->get_extents();
        assert(e.zbegin == 1); assert(e.zend == zsize+1);
        assert(e.ybegin == 1); assert(e.yend == ysize+1);
        assert(e.xbegin == 1); assert(e.xend == xsize+1);
        for (auto zl=e.zbegin; zl<e.zend; ++zl)
            for (auto yl=e.ybegin; yl<e.yend; ++yl)
                for (auto xl=e.xbegin; xl<e.xend; ++xl){
                    auto ndx3d = xl+(yl+zl*ydim)*xdim;
                    sf->at(zl,yl,xl) = static_cast<float>(ndx3d);
                }
        assert(sf->at(1,2,3) == 33);
        assert(sf->at(2,1,3) == 48);
    }

    return 0;
}
