#include "Boundary.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

Boundary::Boundary(const Domain *domain, const LBModel *lbmodel, float solidDensity,
                   BoundaryType type, BoundaryVelocity velocity):
    _domain(domain), _lbmodel(lbmodel), _solidDensity(solidDensity), _type(type), _velocity(velocity){
    std::tie(_xdim, _ydim, _zdim) = domain->getDimensions();
    _kdim = lbmodel->getNumberOfDirections();
}

Boundary::~Boundary(){}

const BoundaryType Boundary::getBoundaryType() const{
    return _type;
}

const BoundaryVelocity Boundary::getBoundaryVelocity() const{
    return _velocity;
}

void Boundary::reset(Lattice &lattice) const{
    applyVelocity(lattice);
    applyDensity(lattice);
}

void Boundary::applyVelocity(Lattice &lattice) const{
    std::vector<std::string> directions = {"east", "west", "north", "south", "up", "down"};
    for (const std::string& dirxn: directions)
        if (_boundaryVelocityIsPrescribed(dirxn))
            _applyVelocityToBoundary(dirxn, lattice);
}

void Boundary::applyDensity(Lattice &lattice) const{
    std::vector<std::string> directions = {"east", "west", "north", "south", "up", "down"};
    for (const std::string& dirxn: directions)
        if (_boundaryTypeIsPrescribed(dirxn))
            if (_type.at(dirxn)=="noslip")
                _applyDensityToBoundary(dirxn, lattice);
}    

void Boundary::applyPeriodicity(Lattice &lattice) const{
    _applyPeriodicityEastWest(lattice);
    _applyPeriodicityNorthSouth(lattice);
}

void Boundary::applyNoslip(Lattice &lattice) const{
    std::vector<std::string> directions = {"east", "west", "north", "south", "up", "down"};
    for (const std::string& dirxn: directions)
        if (_boundaryTypeIsPrescribed(dirxn))
            if (_type.at(dirxn)=="noslip")
                _applyNoslipToBoundary(dirxn, lattice);
}

bool Boundary::_boundaryVelocityIsPrescribed(const std::string direction) const{
    if (_velocity.find(direction)!=_velocity.end())
        return true;
    else
        return false;
}

bool Boundary::_boundaryTypeIsPrescribed(const std::string direction) const{
    if (_type.find(direction)!=_type.end())
        return true;
    else
        return false;
}

const std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>
Boundary::_getBoundaryExtent(const std::string direction) const{
    if (direction=="up"){
        const uint32_t zmin = _zdim, zmax = _zdim;
        const uint32_t ymin = 1, ymax = _ydim;
        const uint32_t xmin = 1, xmax = _xdim;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction=="down"){
        const uint32_t zmin = 1, zmax = 1;
        const uint32_t ymin = 1, ymax = _ydim;
        const uint32_t xmin = 1, xmax = _xdim;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction=="north"){
        const uint32_t zmin = 1, zmax = _zdim;
        const uint32_t ymin = _ydim, ymax = _ydim;
        const uint32_t xmin = 1, xmax = _xdim;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction=="south"){
        const uint32_t zmin = 1, zmax = _zdim;
        const uint32_t ymin = 1, ymax = 1;
        const uint32_t xmin = 1, xmax = _xdim;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction=="east"){
        const uint32_t zmin = 1, zmax = _zdim;
        const uint32_t ymin = 1, ymax = _ydim;
        const uint32_t xmin = _xdim, xmax = _xdim;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    } else if (direction=="west"){
        const uint32_t zmin = 1, zmax = _zdim;
        const uint32_t ymin = 1, ymax = _ydim;
        const uint32_t xmin = 1, xmax = 1;
        return std::tie(xmin, xmax, ymin, ymax, zmin, zmax);
    }
    else{
        throw std::logic_error("Boundary::_getBoundaryExtent: unknow direction "+direction);
    }
}

void Boundary::_applyVelocityToBoundary(const std::string direction, Lattice &lattice) const{
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin, ymax, zmin, zmax) = _getBoundaryExtent(direction);
    for (auto zl=zmin; zl<zmax+1; ++zl)
        for (auto yl=ymin; yl<ymax+1; ++yl)
            for (auto xl=xmin; xl<xmax+1; ++xl)
                for (auto i=0; i<3; ++i)
                    lattice.u->at(zl,yl,xl,i) = _velocity.at(direction)[i];
}

void Boundary::_applyDensityToBoundary(const std::string direction, Lattice &lattice) const{
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin, ymax, zmin, zmax) = _getBoundaryExtent(direction);
    for (auto zl=zmin; zl<zmax+1; ++zl)
        for (auto yl=ymin; yl<ymax+1; ++yl)
            for (auto xl=xmin; xl<xmax+1; ++xl)
                lattice.rho->at(zl,yl,xl) = _solidDensity;
}

void Boundary::_applyNoslipToBoundary(const std::string direction, Lattice &lattice) const{
    const auto c = _lbmodel->getLatticeVelocities();
    const auto w = _lbmodel->getDirectionalWeights();
    const auto reverse = _lbmodel->getReverse();
    const auto cs2inv = 1.0f/_lbmodel->getSpeedOfSoundSquared();
    array4f *n = lattice.n;
    array3f *rho = lattice.rho;
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin, ymax, zmin, zmax) = _getBoundaryExtent(direction);
    auto ub = _velocity.at(direction);
    for (auto zl=zmin; zl<zmax+1; ++zl){
        for (auto yl=ymin; yl<ymax+1; ++yl){
            for (auto xl=xmin; xl<xmax+1; ++xl){
                for (auto kp=1; kp<_kdim; ++kp){ // kp=0 => current node
                    auto kp3 = kp*3;
                    auto nx = xl + c[kp3+0]; // neighbor co-ordinates
                    auto ny = yl + c[kp3+1];
                    auto nz = zl + c[kp3+2];
                    auto rhonbr = rho->at(nz,ny,nx);
                    auto k = reverse[kp]; auto k3 = k*3;
                    auto cu = c[k3+0]*ub[0] + c[k3+1]*ub[1] + c[k3+2]*ub[2];
                    n->at(nz,ny,nx,kp) =
                        n->at(zl,yl,xl,k) - 2.0f*w[k]*rhonbr*cu*cs2inv;
                }
            }
        }
    }               
}

void Boundary::_applyPeriodicityEastWest(Lattice &lattice) const{
    if (_boundaryTypeIsPrescribed("east") && _boundaryTypeIsPrescribed("west")){
        if ((_type.at("east")=="periodic") && (_type.at("west")=="periodic")){
            array4f *n = lattice.n;
            const auto c = _lbmodel->getLatticeVelocities();
            for (auto zl=1; zl<_zdim+1; ++zl){
                for (auto yl=1; yl<_ydim+1; ++yl){
                    for (auto k=0; k<_kdim; ++k){
                        auto k3 = 3*k;
                        // east->west (x-positive components)
                        if (c[k3]>0) n->at(zl,yl,1,k) = n->at(zl,yl,_xdim+1,k);
                        // west->east (x-negative components)
                        if (c[k3]<0) n->at(zl,yl,_xdim,k) = n->at(zl,yl,0,k);
                    }
                }
            }
        }
    }
}

void Boundary::_applyPeriodicityNorthSouth(Lattice &lattice) const{
    if (_boundaryTypeIsPrescribed("north") && _boundaryTypeIsPrescribed("south")){
        if ((_type.at("north")=="periodic") && (_type.at("south")=="periodic")){
            array4f *n = lattice.n;
            const auto c = _lbmodel->getLatticeVelocities();
            for (auto zl=1; zl<_zdim+1; ++zl){
                for (auto xl=1; xl<_xdim+1; ++xl){
                    for (auto k=0; k<_kdim; ++k){
                        auto k3 = 3*k;
                        // north->south (y-positive components)
                        if (c[k3+1]>0) n->at(zl,1,xl,k) = n->at(zl,_ydim+1,xl,k);
                        // south->north (y-negative components)
                        if (c[k3+1]<0) n->at(zl,_ydim,xl,k) = n->at(zl,0,xl,k);
                    }
                }
            }
        }
    }
}
