#include "Boundary.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <chrono>
#include "tbb/tbb.h"

float Boundary::_time_noslip = 0.0f;
float Boundary::_time_periodicity = 0.0f;
float Boundary::_time_reset = 0.0f;

Boundary::Boundary(const LBModel *lbmodel, const Domain *domain, float solid_density,
                   BoundaryType type, BoundaryVelocity velocity):
    _lbmodel(lbmodel), _domain(domain), _solid_density(solid_density), _type(type), _velocity(velocity){
    std::tie(_xdim, _ydim, _zdim) = domain->get_dimensions();
    _kdim = lbmodel->get_num_directions();
}

Boundary::~Boundary(){}

const BoundaryType Boundary::get_boundary_type() const{
    return _type;
}

const BoundaryVelocity Boundary::get_boundary_velocity() const{
    return _velocity;
}

void Boundary::reset(Lattice &lattice) const{
    auto start = std::chrono::system_clock::now();
    
    apply_velocity(lattice);
    apply_density(lattice);

    std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
    Boundary::_time_reset += elapsed.count();
}

void Boundary::apply_velocity(Lattice &lattice) const{
    std::vector<std::string> directions = {"east", "west", "north", "south", "up", "down"};
    for (const std::string& dirxn: directions)
        if (_boundary_velocity_is_prescribed(dirxn))
            _apply_velocity_to_boundary(dirxn, lattice);
}

void Boundary::apply_density(Lattice &lattice) const{
    std::vector<std::string> directions = {"east", "west", "north", "south", "up", "down"};
    for (const std::string& dirxn: directions)
        if (_boundary_type_is_prescribed(dirxn))
            if (_type.at(dirxn)=="noslip")
                _apply_density_to_boundary(dirxn, lattice);
}    

void Boundary::apply_periodicity(Lattice &lattice) const{
    auto start = std::chrono::system_clock::now();

    _apply_periodicity_east_west(lattice);
    _apply_periodicity_north_south(lattice);

    std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
    Boundary::_time_periodicity += elapsed.count();
}

void Boundary::apply_noslip(Lattice &lattice) const{
    auto start = std::chrono::system_clock::now();

    std::vector<std::string> directions = {"east", "west", "north", "south", "up", "down"};
    for (const std::string& dirxn: directions)
        if (_boundary_type_is_prescribed(dirxn))
            if (_type.at(dirxn)=="noslip")
                _apply_noslip_to_boundary(dirxn, lattice);

    std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
    Boundary::_time_noslip += elapsed.count();
}

bool Boundary::_boundary_velocity_is_prescribed(const std::string direction) const{
    if (_velocity.find(direction)!=_velocity.end())
        return true;
    else
        return false;
}

bool Boundary::_boundary_type_is_prescribed(const std::string direction) const{
    if (_type.find(direction)!=_type.end())
        return true;
    else
        return false;
}

const std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>
Boundary::_get_boundary_extent(const std::string direction) const{
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
        throw std::logic_error("Boundary::_get_boundary_extent: unknow direction "+direction);
    }
}

void Boundary::_apply_velocity_to_boundary(const std::string direction, Lattice &lattice) const{
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin, ymax, zmin, zmax) = _get_boundary_extent(direction);

    tbb::parallel_for
        (tbb::blocked_range3d<uint32_t> (zmin, zmax+1, ymin, ymax+1, xmin, xmax+1),
         [this, &lattice, direction]
         (const tbb::blocked_range3d<uint32_t> &r){
            // lambda body - start
            for (auto zl=r.pages().begin(); zl<r.pages().end(); ++zl){
                for (auto yl=r.rows().begin(); yl<r.rows().end(); ++yl){
                    for (auto xl=r.cols().begin(); xl<r.cols().end(); ++xl){
                        for (auto i=0; i<3; ++i)
                            lattice.u->at(zl,yl,xl,i) = _velocity.at(direction)[i];
                    }
                }
            }
            // lambda body - end
        });
                        
    // for (auto zl=zmin; zl<zmax+1; ++zl)
    //     for (auto yl=ymin; yl<ymax+1; ++yl)
    //         for (auto xl=xmin; xl<xmax+1; ++xl)
    //             for (auto i=0; i<3; ++i)
    //                 lattice.u->at(zl,yl,xl,i) = _velocity.at(direction)[i];

}

void Boundary::_apply_density_to_boundary(const std::string direction, Lattice &lattice) const{
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin, ymax, zmin, zmax) = _get_boundary_extent(direction);

    tbb::parallel_for
        (tbb::blocked_range3d<uint32_t> (zmin, zmax+1, ymin, ymax+1, xmin, xmax+1),
         [this, &lattice]
         (const tbb::blocked_range3d<uint32_t> &r){
            // lambda body - start
            for (auto zl=r.pages().begin(); zl<r.pages().end(); ++zl){
                for (auto yl=r.rows().begin(); yl<r.rows().end(); ++yl){
                    for (auto xl=r.cols().begin(); xl<r.cols().end(); ++xl)
                        lattice.rho->at(zl,yl,xl) = _solid_density;
                }
            }
            // lambda body - end
        });
    
    // for (auto zl=zmin; zl<zmax+1; ++zl)
    //     for (auto yl=ymin; yl<ymax+1; ++yl)
    //         for (auto xl=xmin; xl<xmax+1; ++xl)
    //             lattice.rho->at(zl,yl,xl) = _solid_density;

}

void Boundary::_apply_noslip_to_boundary(const std::string direction, Lattice &lattice) const{
    const auto c = _lbmodel->get_lattice_velocities();
    const auto w = _lbmodel->get_directional_weights();
    const auto reverse = _lbmodel->get_reverse();
    const auto cs2inv = 1.0f/_lbmodel->get_speed_of_sound_squared();
    // array4f *n = lattice.n;
    // const array3f *rho = lattice.rho;
    Field::Field<float, 1> *n = lattice.n;
    const Field::Field<float, 1> *rho = lattice.rho;
    size_t xmin, xmax, ymin, ymax, zmin, zmax;
    std::tie(xmin, xmax, ymin, ymax, zmin, zmax) = _get_boundary_extent(direction);
    const auto ub = _velocity.at(direction);
    tbb::parallel_for
        (tbb::blocked_range3d<uint32_t> (zmin, zmax+1, ymin, ymax+1, xmin, xmax+1),
         [this, &c, &w, &reverse, cs2inv, &ub, n, rho]
         (const tbb::blocked_range3d<uint32_t> &r){
            // lambda body - start
            for (auto zl=r.pages().begin(); zl<r.pages().end(); ++zl){
                for (auto yl=r.rows().begin(); yl<r.rows().end(); ++yl){
                    for (auto xl=r.cols().begin(); xl<r.cols().end(); ++xl){
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
            // lambda body -end
        });
}

void Boundary::_apply_periodicity_east_west(Lattice &lattice) const{
    if (_boundary_type_is_prescribed("east") && _boundary_type_is_prescribed("west")){
        if ((_type.at("east")=="periodic") && (_type.at("west")=="periodic")){
            // array4f *n = lattice.n;
            Field::Field<float, 1> *n = lattice.n;
            const auto c = _lbmodel->get_lattice_velocities();
            for (auto zl=1; zl<_zdim+1; ++zl){
                for (auto yl=1; yl<_ydim+1; ++yl){
                    for (auto k=1; k<_kdim; ++k){ // k=0 => rest particle
                        auto ck = &c[k*3];
                        // east->west (x-positive components)
                        if (ck[0]>0) n->at(zl,yl,1,k) = n->at(zl,yl,_xdim+1,k);
                        // west->east (x-negative components)
                        if (ck[0]<0) n->at(zl,yl,_xdim,k) = n->at(zl,yl,0,k);
                    }
                }
            }
        }
    }
}

void Boundary::_apply_periodicity_north_south(Lattice &lattice) const{
    if (_boundary_type_is_prescribed("north") && _boundary_type_is_prescribed("south")){
        if ((_type.at("north")=="periodic") && (_type.at("south")=="periodic")){
            // array4f *n = lattice.n;
            Field::Field<float, 1> *n = lattice.n;
            const auto c = _lbmodel->get_lattice_velocities();
            for (auto zl=1; zl<_zdim+1; ++zl){
                for (auto xl=1; xl<_xdim+1; ++xl){
                    for (auto k=1; k<_kdim; ++k){ // k=0 => rest particle
                        auto ck = &c[k*3];
                        // north->south (y-positive components)
                        if (ck[1]>0) n->at(zl,1,xl,k) = n->at(zl,_ydim+1,xl,k);
                        // south->north (y-negative components)
                        if (ck[1]<0) n->at(zl,_ydim,xl,k) = n->at(zl,0,xl,k);
                    }
                }
            }
        }
    }
}

float Boundary::get_time_noslip() const{
    return _time_noslip;
}

float Boundary::get_time_periodicity() const{
    return _time_periodicity;
}

float Boundary::get_time_reset() const{
    return _time_reset;
}

float Boundary::get_total_time() const{
    return _time_noslip + _time_periodicity + _time_reset;
}
