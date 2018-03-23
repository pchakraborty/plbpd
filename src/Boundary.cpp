#include "Boundary.hpp"

Boundary::Boundary(BoundaryType type, BoundaryVelocity velocity):
    _type(type), _velocity(velocity){}

Boundary::~Boundary(){}

BoundaryType Boundary::getBoundaryType() const{
    return _type;
}

BoundaryVelocity Boundary::getBoundaryVelocity() const{
    return _velocity;
}

void Boundary::apply(Lattice &lattice) const{
}
