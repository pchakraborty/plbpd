#ifndef LB_DYNAMICS_HPP
#define LB_DYNAMICS_HPP

#include <memory>
#include "LBModel.hpp"
#include "Domain.hpp"
#include "CollisionSRT.hpp"
#include "Streaming.hpp"
#include "SimData.hpp"

class LBDynamics final{

private:

    std::shared_ptr<CollisionSRT> _collide;
    std::shared_ptr<Streaming> _stream;

public:

    LBDynamics() = delete;
    LBDynamics(const LBModel* lbmodel, const Domain* domain);
    LBDynamics(const LBModel* lbmodel, const Domain* domain, bool reference);
    LBDynamics(LBDynamics&) = delete;
    LBDynamics& operator=(LBDynamics &) = delete;
    ~LBDynamics();

    void collide(SimData& simdata);
    void stream(SimData& simdata);
    float get_time_collide();
    float get_time_stream();
    float get_total_time();

};

#endif
