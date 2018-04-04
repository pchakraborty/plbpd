#include "Boundary.hpp"

Boundary::Boundary(const Domain &domain, BoundaryType type, BoundaryVelocity velocity):
    _type(type), _velocity(velocity){
    std::tie(_xdim, _ydim, _zdim) = domain.getDimensions();
}

Boundary::~Boundary(){}

const BoundaryType Boundary::getBoundaryType() const{
    return _type;
}

const BoundaryVelocity Boundary::getBoundaryVelocity() const{
    return _velocity;
}

void Boundary::apply(Lattice &lattice) const{
    // 6 bounding faces: E, W, N, S, U, D
    
    if (_velocity.find("east")!=_velocity.end())
        _applyEastBoundary(lattice);
    if (_velocity.find("west")!=_velocity.end())
        _applyWestBoundary(lattice);
    if (_velocity.find("north")!=_velocity.end())
        _applyNorthBoundary(lattice);
    if (_velocity.find("south")!=_velocity.end())
        _applySouthBoundary(lattice);
    if (_velocity.find("up")!=_velocity.end())
        _applyUpBoundary(lattice);
    if (_velocity.find("down")!=_velocity.end())
        _applyDownBoundary(lattice);
}

void Boundary::_applyEastBoundary(Lattice &lattice) const{
    auto xl = _xdim;
    for (auto yl=1; yl<_ydim+1; yl++)
        for (auto zl=1; zl<_zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("east")[i];
}

void Boundary::_applyWestBoundary(Lattice &lattice) const{
    auto xl = 1;
    for (auto yl=1; yl<_ydim+1; yl++)
        for (auto zl=1; zl<_zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("west")[i];
}

void Boundary::_applyNorthBoundary(Lattice &lattice) const{
    auto yl = _ydim;
    for (auto xl=1; xl<_xdim+1; xl++)
        for (auto zl=1; zl<_zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("north")[i];
}

void Boundary::_applySouthBoundary(Lattice &lattice) const{
    auto yl = 1;
    for (auto xl=1; xl<_xdim+1; xl++)
        for (auto zl=1; zl<_zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("south")[i];
}

void Boundary::_applyUpBoundary(Lattice &lattice) const{
    auto zl = _zdim;
    for (auto xl=1; xl<_xdim+1; xl++)
        for (auto yl=1; yl<_ydim+1; yl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("up")[i];
}

void Boundary::_applyDownBoundary(Lattice &lattice) const{
    auto zl = 1;
    for (auto xl=1; xl<_xdim+1; xl++)
        for (auto yl=1; yl<_ydim+1; yl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("down")[i];
}
