#ifndef SRC_SCALARFIELD_HPP_
#define SRC_SCALARFIELD_HPP_

#include "VectorField.hpp"

namespace field {

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

}  // namespace field

#endif  // SRC_SCALARFIELD_HPP_
