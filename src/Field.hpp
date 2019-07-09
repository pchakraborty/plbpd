#ifndef SRC_FIELD_HPP_
#define SRC_FIELD_HPP_

#include <tuple>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <array>

namespace Field {

struct FieldExtents {
    uint32_t zbegin, zend;
    uint32_t ybegin, yend;
    uint32_t xbegin, xend;
};

template <typename T, uint32_t num_buffer_layers>
class VectorField {
 private:
    uint32_t _zlen, _ylen, _xlen;  // field dimensions
    uint32_t _vlen;  // vector length at each field location
    std::vector<T> _arrdata;
    uint32_t _zfactor, _yfactor, _xfactor;  // factors for index calculation
    FieldExtents _e;

    void _check_non_zero_arg(const uint32_t arg) {
        if (arg == 0)
            throw std::invalid_argument("argument is zero");
    }

    void _set_lengths(const uint32_t zlen,
                      const uint32_t ylen,
                      const uint32_t xlen,
                      const uint32_t vlen) {
        _zlen = zlen + 2*num_buffer_layers;
        _ylen = ylen + 2*num_buffer_layers;
        _xlen = xlen + 2*num_buffer_layers;
        _vlen = vlen;
    }

    void _set_factors() {
        _zfactor = _vlen*_xlen*_ylen;
        _yfactor = _vlen*_xlen;
        _xfactor = _vlen;
    }

    void _set_extents() {
        _e.zbegin = 0 + num_buffer_layers;
        _e.zend = _zlen - num_buffer_layers;
        _e.ybegin = 0 + num_buffer_layers;
        _e.yend = _ylen - num_buffer_layers;
        _e.xbegin = 0 + num_buffer_layers;
        _e.xend = _xlen - num_buffer_layers;
    }

 public:
    // Initialize a field of vectors (of length vlen)
    VectorField(const uint32_t zlen,
                const uint32_t ylen,
                const uint32_t xlen,
                const uint32_t vlen,
                const T& value) {
        for (auto arg : std::vector<uint32_t> {zlen, ylen, xlen, vlen})
            _check_non_zero_arg(arg);
        _set_lengths(zlen, ylen, xlen, vlen);
        _arrdata.resize(_zlen*_ylen*_xlen*_vlen, value);
        _set_factors();
        _set_extents();
    }

    ~VectorField() {}

    VectorField(VectorField&) = delete;

    VectorField& operator=(VectorField&) = delete;

    inline std::tuple<uint32_t, uint32_t, uint32_t> get_dimensions() const {
        return std::make_tuple(_zlen, _ylen, _xlen);
    }

    inline uint32_t get_vector_length() const {
        return _vlen;
    }

    inline const FieldExtents& get_extents() const {
        return _e;
    }

    inline uint32_t sub2ind(uint32_t z, uint32_t y, uint32_t x, uint32_t v) {
        return z*_zfactor + y*_yfactor + x*_xfactor + v;
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v) {
        return _arrdata[sub2ind(z, y, x, v)];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const {
        return _arrdata[sub2ind(z, y, x, v)];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const {
        return &_arrdata[sub2ind(z, y, x, v)];
    }

    inline T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v) {
        return &_arrdata[sub2ind(z, y, x, v)];
    }

    inline std::tuple<uint32_t, uint32_t, uint32_t>
    get_neighbor(uint32_t z, uint32_t y, uint32_t x,
                 const std::array<int32_t, 3> &ck) const {
        return std::make_tuple(z+ck[2], y+ck[1], x+ck[0]);
    }
};  // class VectorField

template <typename T, uint32_t num_buffer_layers>
class ScalarField : public VectorField<T, num_buffer_layers> {
 public:
    ScalarField(uint32_t zlen, uint32_t ylen, uint32_t xlen, const T& value)
        : VectorField<T, num_buffer_layers>(zlen, ylen, xlen, 1, value) {}

    ~ScalarField() {}

    ScalarField(ScalarField&) = delete;

    ScalarField& operator=(ScalarField&) = delete;

    inline uint32_t sub2ind(uint32_t z, uint32_t y, uint32_t x) {
        return VectorField<T, num_buffer_layers>::sub2ind(z, y, x, 0);
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x) {
        return VectorField<T, num_buffer_layers>::at(z, y, x, 0);
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x) const {
        return VectorField<T, num_buffer_layers>::at(z, y, x, 0);
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x) const {
        return VectorField<T, num_buffer_layers>::get(z, y, x, 0);
    }

    inline T* get(uint32_t z, uint32_t y, uint32_t x) {
        return VectorField<T, num_buffer_layers>::get(z, y, x, 0);
    }
};  // class SclalarField

}  // namespace Field

#endif  // SRC_FIELD_HPP_
