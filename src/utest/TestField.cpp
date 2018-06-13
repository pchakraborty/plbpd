#include <iostream>
#include "../Field.hpp"
#undef NDEBUG
#include <cassert>
#include <vector>
#include <memory>

void assert_extents(uint32_t zsize, uint32_t ysize, uint32_t xsize,
                    uint32_t nlayers, const Field::FieldExtents& e){
    assert(e.zbegin == 0+nlayers); assert(e.zend == zsize+nlayers);
    assert(e.ybegin == 0+nlayers); assert(e.yend == ysize+nlayers);
    assert(e.xbegin == 0+nlayers); assert(e.xend == xsize+nlayers);
}

void test_scalar_field_with_0_buffer_layers(){
    auto zsize = 4;
    auto ysize = 5;
    auto xsize = 7;

    // Scalar field
    Field::Field<float,0> sf(zsize, ysize, xsize, 0);
    assert(sf.get_dimensions() == std::tie(zsize, ysize, xsize));
    assert(sf.get_vector_length() == 1);

    // Set data using iterator
    auto ctr = 0;
    for (auto it=sf.begin(); it!=sf.end(); ++it)
        *it = ctr++;

    // Get data using at(z,y,x)
    Field::FieldExtents e = sf.get_extents();
    assert_extents(zsize, ysize, xsize, 0, e);
    for (auto zl=e.zbegin; zl<e.zend; ++zl)
        for (auto yl=e.ybegin; yl<e.yend; ++yl)
            for (auto xl=e.xbegin; xl<e.xend; ++xl)
                assert(sf.at(zl,yl,xl) == sf.get_linear_index(zl,yl,xl));
}

void test_vector_field_with_0_buffer_layers(){
    // Vector field with 0 buffer layers

    const auto zsize = 4;
    const auto ysize = 5;
    const auto xsize = 7;
    const auto vsize = 2;

    Field::Field<float,0> vf(zsize, ysize, xsize, vsize, 0.0f);
    assert(vf.get_dimensions() == std::tie(zsize, ysize, xsize));
    assert(vf.get_vector_length() == vsize);
    // Set data using iterator
    auto ctr = 0;
    for (auto it=vf.begin(); it!=vf.end(); ++it)
        *it = ctr++;

    // Get data using at(z,y,x)
    Field::FieldExtents e = vf.get_extents();
    assert_extents(zsize, ysize, xsize, 0, e);
    auto vlen = vf.get_vector_length();
    for (auto zl=e.zbegin; zl<e.zend; ++zl)
        for (auto yl=e.ybegin; yl<e.yend; ++yl)
            for (auto xl=e.xbegin; xl<e.xend; ++xl){
                auto fndx = vf.get_linear_index(zl,yl,xl,0); // field index
                for (auto kl=0; kl<vsize; ++kl)
                    assert(vf.at(zl,yl,xl,kl) == kl+fndx);
            }
}

void test_scalar_field_with_1_buffer_layers(){
    // Vector field with 1 buffer layers
    const auto zsize=2;
    const auto ysize=2;
    const auto xsize=3;

    // Scalar field with one buffer layer
    auto sf = std::make_unique<Field::Field<float,1> >(zsize, ysize, xsize, 0.0f);
    size_t zdim, ydim, xdim;
    std::tie(zdim, ydim, xdim) = sf->get_dimensions();
    assert(zdim == zsize+2); assert(ydim == ysize+2); assert(xdim == xsize+2);
    Field::FieldExtents e = sf->get_extents();
    assert_extents(zsize, ysize, xsize, 1, e);
    for (auto zl=e.zbegin; zl<e.zend; ++zl)
        for (auto yl=e.ybegin; yl<e.yend; ++yl)
            for (auto xl=e.xbegin; xl<e.xend; ++xl)
                sf->at(zl,yl,xl) = static_cast<float>(sf->get_linear_index(zl,yl,xl));
    assert(sf->at(1,2,3) == 33);
    assert(sf->at(2,1,3) == 48);
}

int main(){

    test_scalar_field_with_0_buffer_layers();
    test_vector_field_with_0_buffer_layers();
    test_scalar_field_with_1_buffer_layers();

    return 0;
}
