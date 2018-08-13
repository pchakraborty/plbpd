#ifndef PLBPD_HPP
#define PLBPD_HPP

#include "Flow.hpp"
#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"
#include "SimData.hpp"
#include "LBDynamics.hpp"
#include "CalcMoments.hpp"

class PLBPD final {
 private:
    const LBModel *_lbmodel;
    const Domain *_domain;
    const Boundary *_boundary;
    std::unique_ptr<LBDynamics> _lbdynamics;
    std::unique_ptr<CalcMoments> _calc_moments;

 public:
    PLBPD(const LBModel *lbmodel,
          const Domain *domain,
          const Boundary *boundary);
    PLBPD(PLBPD&) = delete;
    PLBPD& operator=(PLBPD&) = delete;
    ~PLBPD();
    void initialize(SimData& simdata) const;
    void march_in_time(size_t num_timesteps, SimData& simdata) const;
    void print_times() const;
};

#endif
