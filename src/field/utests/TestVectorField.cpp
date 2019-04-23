#include <cassert>
#include <tuple>

#include "../VectorField.hpp"

template<typename T, uint32_t n_halo_layers>
void assert_extents(field::VectorField<T, n_halo_layers>& vf) {
    field::FieldExtent e = vf.get_field_extent();
    uint32_t zsize, ysize, xsize;
    std::tie(zsize, ysize, xsize) = vf.get_field_size();
    assert(e.zbegin == 0 + n_halo_layers);
    assert(e.zend == zsize + n_halo_layers);
    assert(e.ybegin == 0 + n_halo_layers);
    assert(e.yend == ysize + n_halo_layers);
    assert(e.xbegin == 0 + n_halo_layers);
    assert(e.xend == xsize + n_halo_layers);
}

void test_vectorfield_with_0_halo_layers() {
    const uint32_t zsize = 4, ysize = 5, xsize = 3; uint32_t vsize = 4;
    const uint32_t n_halo_layers = 0;
    
    field::VectorField<float, n_halo_layers> vf(zsize, ysize, xsize, vsize);
    
    assert(std::make_tuple(zsize, ysize, xsize) == vf.get_field_size());
    assert(vsize == vf.get_vector_length());
    assert_extents<float, n_halo_layers>(vf);    
}

int main() {
    test_vectorfield_with_0_halo_layers();
    
    return 0;
}
