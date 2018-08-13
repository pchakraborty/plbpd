#include "PLBPD.hpp"

#include <iostream>
#include "tbb/tbb.h"

PLBPD::PLBPD(const LBModel *lbmodel,
             const Domain *domain,
             const Boundary *boundary)
    : _lbmodel(lbmodel), _domain(domain), _boundary(boundary) {
    // Dynamics
    _lbdynamics = std::make_unique<LBDynamics>(_lbmodel, _domain);
    // Moment calculator
    _calc_moments = std::make_unique<CalcMoments>(_lbmodel);
}

void PLBPD::initialize(SimData& simdata) const {
    _domain->initialize(simdata);
    _boundary->reset(simdata);
}

void PLBPD::march_in_time(size_t num_timesteps, SimData& simdata) const {
    tbb::tick_count start = tbb::tick_count::now();
    for (auto istep = 0; istep < num_timesteps; ++istep) {
        _lbdynamics->collide(simdata);
        _lbdynamics->stream(simdata, "push");
        _boundary->apply_noslip(simdata);
        _boundary->apply_periodicity(simdata);
        (*_calc_moments)(simdata);
        _boundary->reset(simdata);
    }
    auto elapsed = (tbb::tick_count::now()-start).seconds();
    std::cout << "Time taken by timeloop: " << elapsed << "s" << std::endl;
}

PLBPD::~PLBPD() {}

void PLBPD::print_times() const {
    std::cout << "LBDynamics: " << _lbdynamics->get_total_time() << "s\n";
    std::cout << "-collide: " << _lbdynamics->get_time_collide() << "s\n";
    std::cout << "-stream: " << _lbdynamics->get_time_stream() << "s\n";
    std::cout << "CalcMoments: " << _calc_moments->get_total_time() << "s\n";
    std::cout << "Boundary: " << _boundary->get_total_time() << "s\n";
    std::cout << "-noslip: " << _boundary->get_time_noslip() << "s\n";
    std::cout << "-periodicity: " << _boundary->get_time_periodicity() << "s\n";
    std::cout << "-reset: " << _boundary->get_time_reset() << "s\n";
}
