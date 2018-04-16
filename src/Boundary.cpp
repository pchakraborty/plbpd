#include "Boundary.hpp"
#include <iostream>

Boundary::Boundary(const Domain *domain, const LBModel *lbmodel, BoundaryType type, BoundaryVelocity velocity):
    _type(type), _velocity(velocity){
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

void Boundary::apply(Lattice &lattice) const{
    // 6 bounding faces: E, W, N, S, U, D
    // if (_velocity.find("east")!=_velocity.end())
    //     _applyEastBoundary(lattice);
    // if (_velocity.find("west")!=_velocity.end())
    //     _applyWestBoundary(lattice);
    // if (_velocity.find("north")!=_velocity.end())
    //     _applyNorthBoundary(lattice);
    // if (_velocity.find("south")!=_velocity.end())
    //     _applySouthBoundary(lattice);
    _applyUpBoundary(lattice);
    _applyDownBoundary(lattice);
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


void Boundary::_applyEastBoundary(Lattice &lattice) const{
    // std::cout<<"east\n";
    // const auto xl = _xdim;
    // for (auto yl=1; yl<_ydim+1; ++yl)
    //     for (auto zl=1; zl<_zdim+1; ++zl)
    //         for (auto i=0; i<3; i++)
    //             lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("east")[i];
}

void Boundary::_applyWestBoundary(Lattice &lattice) const{
    // std::cout<<"west\n";
    // const auto xl = 1;
    // for (auto yl=1; yl<_ydim+1; ++yl)
    //     for (auto zl=1; zl<_zdim+1; ++zl)
    //         for (auto i=0; i<3; i++)
    //             lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("west")[i];
}

void Boundary::_applyNorthBoundary(Lattice &lattice) const{
    // std::cout<<"north\n";
    // const auto yl = _ydim;
    // for (auto xl=1; xl<_xdim+1; ++xl)
    //     for (auto zl=1; zl<_zdim+1; ++zl)
    //         for (auto i=0; i<3; i++)
    //             lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("north")[i];
}

void Boundary::_applySouthBoundary(Lattice &lattice) const{
    // std::cout<<"south\n";
    // const auto yl = 1;
    // for (auto xl=1; xl<_xdim+1; ++xl)
    //     for (auto zl=1; zl<_zdim+1; ++zl)
    //         for (auto i=0; i<3; i++)
    //             lattice.u[i+(zl+(yl+xl*_ydim)*_zdim)*3] = _velocity.at("south")[i];
}

void Boundary::_applyUpBoundary(Lattice &lattice) const{
    const auto zl = _zdim;

    // Boundary velocity
    if (_boundaryVelocityIsPrescribed("up"))
        for (auto yl=1; yl<_ydim+1; ++yl)
            for (auto xl=1; xl<_xdim+1; ++xl)
                for (auto i=0; i<3; ++i)
                    lattice.u->at(zl,yl,xl,i) = _velocity.at("up")[i];

    // // Reflection from solid boundary
    // const auto c = _lbmodel->getLatticeVelocities();
    // if (_boundaryTypeIsPrescribed("down")){ // down bdry type is prescribed
    //     if (_type.at("down")=="noslip"){
    //     for (auto yl=1; yl<_ydim+1; ++yl){
    //         auto tmp = (yl+zl*(_ydim+2))*(_xdim+2);
    //         for (auto xl=1; xl<_xdim+1; ++xl){
    //             auto zyx = xl+tmp;
    //             for (auto k=1; k<_kdim; ++k){ // k=0 => current node
    //                 auto k3 = k*3;
    //                 std::array<int, 3> ck = {c[k3], c[k3+1], c[k3+2]};
    //                 auto knbr = ;
    //             }
    //         }
    //     }
    // }
    
}

void Boundary::_applyDownBoundary(Lattice &lattice) const{
    const auto zl = 1;

    // Boundary velocity
    if (_boundaryVelocityIsPrescribed("down"))
        for (auto yl=1; yl<_ydim+1; ++yl)
            for (auto xl=1; xl<_xdim+1; ++xl)
                for (auto i=0; i<3; ++i)
                    lattice.u->at(zl,yl,xl,i) = _velocity.at("down")[i];
    
    // // Reflection from solid boundary
    // const auto c = _lbmodel->getLatticeVelocities();
    // if (_boundaryTypeIsPrescribed("down")){ // down bdry type is prescribed
    //     if (_type.at("down")=="noslip"){
    //     for (auto yl=1; yl<_ydim+1; ++yl){
    //         auto tmp = (yl+zl*(_ydim+2))*(_xdim+2);
    //         for (auto xl=1; xl<_xdim+1; ++xl){
    //             auto zyx = xl+tmp;
    //             for (auto k=1; k<_kdim; ++k){ // k=0 => current node
    //                 auto k3 = k*3;
    //                 std::array<int, 3> ck = {c[k3], c[k3+1], c[k3+2]};
    //                 auto knbr = ;
    //             }
    //         }
    //     }
    // }

}
