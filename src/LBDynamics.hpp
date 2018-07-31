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
    void _init(
        const LBModel* lbmodel,
        const Domain* domain,
        std::string stream_type,
        bool reference
    );

public:

    LBDynamics() = delete;
    LBDynamics(const LBModel* lbmodel, const Domain* domain);
    LBDynamics(const LBModel* lbmodel, const Domain* domain, bool reference);
    LBDynamics(const LBModel* lbmodel, const Domain* domain, std::string stream_type);
    LBDynamics(
        const LBModel* lbmodel,
        const Domain* domain,
        std::string stream_type,
        bool reference
    );
    LBDynamics(LBDynamics&) = delete;
    LBDynamics& operator=(LBDynamics &) = delete;
    ~LBDynamics();

    void collide(SimData& simdata) const;
    void stream(SimData& simdata) const;
    float get_time_collide() const;
    float get_time_stream() const;
    float get_total_time() const;

};

#endif
