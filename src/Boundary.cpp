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

void Boundary::apply(const Domain &domain, Lattice &lattice){
    // 6 bounding faces: E, W, N, S, U, D
    
    if (_velocity.find("east")!=_velocity.end())
        _applyEastBoundary(domain, lattice);
    if (_velocity.find("west")!=_velocity.end())
        _applyWestBoundary(domain, lattice);
    if (_velocity.find("north")!=_velocity.end())
        _applyNorthBoundary(domain, lattice);
    if (_velocity.find("south")!=_velocity.end())
        _applySouthBoundary(domain, lattice);
    if (_velocity.find("up")!=_velocity.end())
        _applyUpBoundary(domain, lattice);
    if (_velocity.find("down")!=_velocity.end())
        _applyDownBoundary(domain, lattice);
}

void Boundary::_applyEastBoundary(const Domain &domain, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = domain.getDomainDimensions();
    auto xl = xdim;
    for (auto yl=1; yl<ydim+1; yl++)
        for (auto zl=1; zl<zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*ydim)*zdim)*3] = _velocity["east"][i];
}

void Boundary::_applyWestBoundary(const Domain &domain, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = domain.getDomainDimensions();
    auto xl = 1;
    for (auto yl=1; yl<ydim+1; yl++)
        for (auto zl=1; zl<zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*ydim)*zdim)*3] = _velocity["west"][i];
}

void Boundary::_applyNorthBoundary(const Domain &domain, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = domain.getDomainDimensions();
    auto yl = ydim;
    for (auto xl=1; xl<xdim+1; xl++)
        for (auto zl=1; zl<zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*ydim)*zdim)*3] = _velocity["north"][i];
}

void Boundary::_applySouthBoundary(const Domain &domain, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = domain.getDomainDimensions();
    auto yl = 1;
    for (auto xl=1; xl<xdim+1; xl++)
        for (auto zl=1; zl<zdim+1; zl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*ydim)*zdim)*3] = _velocity["south"][i];
}

void Boundary::_applyUpBoundary(const Domain &domain, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = domain.getDomainDimensions();
    auto zl = zdim;
    for (auto xl=1; xl<xdim+1; xl++)
        for (auto yl=1; yl<ydim+1; yl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*ydim)*zdim)*3] = _velocity["up"][i];
}

void Boundary::_applyDownBoundary(const Domain &domain, Lattice &lattice){
    size_t xdim, ydim, zdim;
    std::tie(xdim, ydim, zdim) = domain.getDomainDimensions();
    auto zl = 1;
    for (auto xl=1; xl<xdim+1; xl++)
        for (auto yl=1; yl<ydim+1; yl++)
            for (auto i=0; i<3; i++)
                lattice.u[i+(zl+(yl+xl*ydim)*zdim)*3] = _velocity["down"][i];
}
