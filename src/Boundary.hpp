#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <string>
#include <unordered_map>
#include "Lattice.hpp"

using BoundaryType = std::unordered_map<std::string, std::string>;
using BoundaryVelocity = std::unordered_map<std::string, std::array<float, 3> >;

class Boundary final{

private:

    BoundaryType _type;
    BoundaryVelocity _velocity;

public:

    Boundary(BoundaryType type, BoundaryVelocity velocity);
    ~Boundary();
    BoundaryType getBoundaryType() const;
    BoundaryVelocity getBoundaryVelocity() const;
    void apply(Lattice &lattice) const;

};

#endif
