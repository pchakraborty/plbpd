#include <cassert>
#include <tuple>
#include <iostream>

#include "../VectorField.hpp"

template<typename T, uint32_t n_halo_layers>
void assert_extents(field::VectorField<T, n_halo_layers>& vf) {
    field::FieldExtent e = vf.get_extent();
    std::array<uint32_t, 3> shape = vf.get_shape();
    assert(e.zbegin == 0 + n_halo_layers);
    assert(e.zend == shape[0] + n_halo_layers);
    assert(e.ybegin == 0 + n_halo_layers);
    assert(e.yend == shape[1] + n_halo_layers);
    assert(e.xbegin == 0 + n_halo_layers);
    assert(e.xend == shape[2] + n_halo_layers);
}

void test_vectorfield_with_0_halo_layers() {
    const std::array<uint32_t, 3>shape{4, 5, 3};
    const uint32_t vlen = 2;
    const uint32_t n_halo_layers = 0;
    
    field::VectorField<float, n_halo_layers> vf(shape[0], shape[1], shape[2], vlen);
    
    assert(shape == vf.get_shape());
    assert(vlen == vf.get_vector_length());
    assert_extents<float, n_halo_layers>(vf);    
}

int main() {
    test_vectorfield_with_0_halo_layers();
    
    return 0;
}
