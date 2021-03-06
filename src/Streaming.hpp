#ifndef SRC_STREAMING_HPP_
#define SRC_STREAMING_HPP_

#include <string>
#include <vector>

#include "SimData.hpp"
#include "LBModel.hpp"

class Streaming final {
 private:
    const std::vector<std::array<int32_t, 3> > &_c;  // directional velocities
    static float _time;

    void _push_ref(SimData &simdata) const;
    void _push_tbb(SimData &simdata) const;
    void _pull_ref(SimData &simdata) const;
    void _pull_tbb(SimData &simdata) const;
    void _stream_ref(SimData &simdata, std::string stream_type) const;
    void _stream_tbb(SimData &simdata, std::string stream_type) const;
    void _push_kernel(
        const size_t zl, const size_t yl, const size_t xl,
        SimData &simdata) const;
 public:
    explicit Streaming(const LBModel *lbmodel);
    Streaming(Streaming&) = delete;
    Streaming& operator=(Streaming&) = delete;
    ~Streaming();
    void operator()(
        SimData &simdata,
        std::string stream_type,
        bool reference = false) const;
    float get_total_time() const;
};

#endif  // SRC_STREAMING_HPP_
