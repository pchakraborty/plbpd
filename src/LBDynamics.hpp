#ifndef SRC_LBDYNAMICS_HPP_
#define SRC_LBDYNAMICS_HPP_

#include <memory>
#include <string>

#include "LBModel.hpp"
#include "Domain.hpp"
#include "CollisionSRT.hpp"
#include "Streaming.hpp"
#include "SimData.hpp"

class LBDynamics final {
 private:
    std::unique_ptr<CollisionSRT> _collide;
    std::unique_ptr<Streaming> _stream;

 public:
    LBDynamics() = delete;
    LBDynamics(const LBModel* lbmodel, const Domain* domain);
    LBDynamics(LBDynamics&) = delete;
    LBDynamics& operator=(LBDynamics &) = delete;
    ~LBDynamics();

    void collide(SimData& simdata, bool reference = false) const;
    void stream(
        SimData& simdata,
        std::string stream_type,
        bool reference = false) const;
    float get_time_collide() const;
    float get_time_stream() const;
    float get_total_time() const;
};

#endif  // SRC_LBDYNAMICS_HPP_
