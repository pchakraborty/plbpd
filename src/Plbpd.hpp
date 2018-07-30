#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include "Flow.hpp"
#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "SimData.hpp"
#include "LBDynamics.hpp"
#include "CalcMoments.hpp"

class Plbpd final{
    
private:

    const LBModel *_lbmodel;
    const Domain *_domain;
    const Boundary *_boundary;
    std::shared_ptr<LBDynamics> _lbdynamics;
    std::shared_ptr<CalcMoments> _calc_moments;
    
public:

    Plbpd(const LBModel *lbmodel, const Domain *domain, const Boundary *boundary);
    Plbpd(Plbpd&) = delete;
    Plbpd& operator=(Plbpd&) = delete;
    ~Plbpd();
    void initialize(SimData& simdata) const;
    void march_in_time(size_t num_timesteps, SimData& simdata) const;
    void print_times() const;
    
};

#endif
