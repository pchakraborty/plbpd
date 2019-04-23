#ifndef SRC_VECTORFIELD_HPP_
#define SRC_VECTORFIELD_HPP_

#include <tuple>
#include <vector>
#include <stdexcept>
#include <memory>

#include "FieldExtent.hpp"

namespace field {

template <typename T, uint32_t n_halo_layers>
class VectorField {

    using uintArray3 = std::array<uint32_t, 3>;
    using intArray3 = std::array<int32_t, 3>;

public:

    // Initialize a field of vectors (of length vsize)
    VectorField(uint32_t zsize, uint32_t ysize, uint32_t xsize, uint32_t vlen) {
        // _sanity_check(std::array<uint32_t, 4>{zsize, ysize, xsize, vlen};
        _check_positive_arg(zsize);
        _check_positive_arg(ysize);
        _check_positive_arg(xsize);
        _check_positive_arg(vlen);
        _set_lengths(zsize, ysize, xsize, vlen);
        _arrdata.resize(_zsize*_ysize*_xsize*_vsize, static_cast<T>(0));
        _set_factors();
        _set_extents();
    }

    ~VectorField() {}

    VectorField(VectorField&) = delete;
    VectorField& operator=(VectorField&) = delete;

    inline std::tuple<uint32_t, uint32_t, uint32_t> get_field_size() const {
        return std::make_tuple(_zsize, _ysize, _xsize);
    }

    inline uint32_t get_vector_length() const {
        return _vsize;
    }

    inline const FieldExtent& get_field_extent() const {
        return _e;
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v) {
        return _arrdata[_sub2ind(z, y, x, v)];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const {
        return _arrdata[_sub2ind(z, y, x, v)];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const {
        return &_arrdata[_sub2ind(z, y, x, v)];
    }

    inline T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v) {
        return &_arrdata[_sub2ind(z, y, x, v)];
    }

    inline intArray3
    get_neighbor(uint32_t z, uint32_t y, uint32_t x,  const intArray3 & ck) const {
        return intArray3{z+ck[2], y+ck[1], x+ck[0]};
    }

private:

    uint32_t _zsize, _ysize, _xsize;  // field size in each direction
    uint32_t _vsize;  // vector length at each field location
    std::vector<T> _arrdata;
    FieldExtent _e;
    uint32_t _zfactor, _yfactor, _xfactor;  // factors for index calculation

    void _check_positive_arg(const uint32_t arg) {
        if (arg < 1)  throw std::invalid_argument("argument is less than 1");
    }

    void _set_lengths(uint32_t zsize, uint32_t ysize, uint32_t xsize, uint32_t vsize) {
        _zsize = zsize + 2*n_halo_layers;
        _ysize = ysize + 2*n_halo_layers;
        _xsize = xsize + 2*n_halo_layers;
        _vsize = vsize;
    }

    void _set_factors() {
        _zfactor = _vsize*_xsize*_ysize;
        _yfactor = _vsize*_xsize;
        _xfactor = _vsize;
    }

    void _set_extents() {
        _e.zbegin = 0 + n_halo_layers;
        _e.zend = _zsize - n_halo_layers;
        _e.ybegin = 0 + n_halo_layers;
        _e.yend = _ysize - n_halo_layers;
        _e.xbegin = 0 + n_halo_layers;
        _e.xend = _xsize - n_halo_layers;
    }

    inline uint32_t _sub2ind(uint32_t z, uint32_t y, uint32_t x, uint32_t v) {
        return z*_zfactor + y*_yfactor + x*_xfactor + v;
    }


};  // class VectorField

}  // namespace field

#endif  // SRC_VECTORFIELD_HPP_
