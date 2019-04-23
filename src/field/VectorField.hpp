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

public:

    // Initialize a field of vectors (of length vlen)
    VectorField(std::array<uint32_t, 3> field_shape, uint32_t vlen) {
        _initialize(field_shape, vlen);
    }

    VectorField(uint32_t shp1, uint32_t shp2, uint32_t shp3, uint32_t vlen) {
        std::array<uint32_t, 3> field_shape{shp1, shp2, shp3};
        _initialize(field_shape, vlen);
    }

    ~VectorField() {}

    VectorField(VectorField&) = delete;
    VectorField& operator=(VectorField&) = delete;

    inline std::array<uint32_t, 3> get_shape() const {
        return std::array<uint32_t, 3>{_shape1, _shape2, _shape3};
    }

    inline uint32_t get_vector_length() const {
        return _vlen;
    }

    inline const FieldExtent& get_extent() const {
        return _e;
    }

    inline T& at(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t v) {
        return _arrdata[_sub2ind(x1, x2, x3, v)];
    }

    inline const T& at(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t v) const {
        return _arrdata[_sub2ind(x1, x2, x3, v)];
    }

    inline const T* get(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t v) const {
        return &_arrdata[_sub2ind(x1, x2, x3, v)];
    }

    inline T* get(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t v) {
        return &_arrdata[_sub2ind(x1, x2, x3, v)];
    }

    inline std::array<int32_t, 3> get_neighbor(
        uint32_t x1, uint32_t x2, uint32_t x3,
        const std::array<int32_t, 3> & ck) const
    {
        return std::array<int32_t, 3>{x1+ck[2], x2+ck[1], x3+ck[0]};
    }

private:

    uint32_t _shape1, _shape2, _shape3;  // field size in each direction
    uint32_t _vlen;  // vector length at each field location
    std::vector<T> _arrdata;
    FieldExtent _e;
    uint32_t _factor1, _factor2, _factor3;  // factors for index calculation

    void _initialize(std::array<uint32_t, 3>& field_shape, uint32_t vlen) {
        _set_lengths(field_shape, vlen);
        _create_internal_array();
        _set_factors();
        _set_extents();
    }        

    void _check_positive_arg(const uint32_t arg) {
        if (arg < 1)  throw std::invalid_argument("argument is less than 1");
    }

    void _set_lengths(std::array<uint32_t, 3>& shape, uint32_t vlen) {
        for (auto& element: shape)
            if (element < 1)
                throw std::invalid_argument("Field shape cannot be < 1");
        _shape1 = shape[0] + 2*n_halo_layers;
        _shape2 = shape[1] + 2*n_halo_layers;
        _shape3 = shape[2] + 2*n_halo_layers;
        if (vlen < 0)
            throw std::invalid_argument("Vector length cannot be < 0");
        _vlen = vlen;
    }

    void _create_internal_array() {
        _arrdata.resize(_shape1*_shape2*_shape3*_vlen, static_cast<T>(0));
    }

    void _set_factors() {
        _factor1 = _vlen*_shape3*_shape2;
        _factor2 = _vlen*_shape3;
        _factor3 = _vlen;
    }

    void _set_extents() {
        _e.zbegin = 0 + n_halo_layers;
        _e.zend = _shape1 - n_halo_layers;
        _e.ybegin = 0 + n_halo_layers;
        _e.yend = _shape2 - n_halo_layers;
        _e.xbegin = 0 + n_halo_layers;
        _e.xend = _shape3 - n_halo_layers;
    }

    inline uint32_t _sub2ind(uint32_t z, uint32_t y, uint32_t x, uint32_t v) {
        return z*_factor1 + y*_factor2 + x*_factor3 + v;
    }


};  // class VectorField

}  // namespace field

#endif  // SRC_VECTORFIELD_HPP_
