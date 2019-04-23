#ifndef SRC_FIELDEXTENTS_HPP_
#define SRC_FIELDEXTENTS_HPP_

namespace field {

struct FieldExtent {
    uint32_t zbegin, zend;
    uint32_t ybegin, yend;
    uint32_t xbegin, xend;
};

}  // namespace field

#endif  // SRC_FIELDEXTENTS_HPP_
