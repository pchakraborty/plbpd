#include "Streaming.hpp"

#include <chrono>
#include <string>
#include <algorithm>
#include <vector>
#include "tbb/tbb.h"

namespace chrono = std::chrono;

float Streaming::_time = 0.0f;

float Streaming::get_total_time() const {
    return _time;
}

bool is_valid(std::string stream_type) {
    std::vector<std::string> types = {"push", "pull"};
    if (std::find(types.begin(), types.end(), stream_type) != types.end())
        return true;
    else
        throw std::logic_error("Invalid stream type: " + stream_type);
}

Streaming::Streaming(const LBModel *lbmodel)
    : _c(lbmodel->get_directional_velocities()) {}

Streaming::~Streaming() {}

void Streaming::operator()(
    SimData &simdata,
    std::string stream_type,
    bool reference /*= false*/) const {
    auto start = chrono::system_clock::now();
    if (reference)
        _stream_ref(simdata, stream_type);
    else
        _stream_tbb(simdata, stream_type);
    chrono::duration<float> elapsed = chrono::system_clock::now()-start;
    Streaming::_time += elapsed.count();
}

void Streaming::_stream_ref(SimData &simdata, std::string stream_type) const {
    if (is_valid(stream_type)) {
        if (stream_type == "push")
            _push_ref(simdata);
        else // stream_type == pull
            _pull_ref(simdata);
    }
}

void Streaming::_stream_tbb(SimData &simdata, std::string stream_type) const {
    if (is_valid(stream_type)) {
        if (stream_type == "push")
            _push_tbb(simdata);
        else // stream_type == pull
            _pull_tbb(simdata);
    }
}

__attribute__((always_inline))
inline void Streaming::_push_kernel(
    size_t zl, size_t yl, size_t xl,
    SimData &simdata) const {
    // Streaming (push) kernel
    auto nlocal = simdata.n->get(zl, yl, xl, 0);
    auto kdim = simdata.n->get_vector_length();
    for (auto k = 0; k < kdim; ++k) {
        auto ck = &_c[k*3];
        size_t nz, ny, nx;  // k-nbr of (zl, yl, xl)
        std::tie(nz, ny, nx) =  simdata.n->get_neighbor(zl, yl, xl, ck);
        simdata.ntmp->at(nz, ny, nx, k) = nlocal[k];
    }
}

void Streaming::_push_ref(SimData &simdata) const {
    const auto e = simdata.n->get_extents();
    for (auto zl = e.zbegin; zl < e.zend; ++zl)
        for (auto yl = e.ybegin; yl < e.yend; ++yl)
            for (auto xl = e.xbegin; xl < e.xend; ++xl)
                _push_kernel(zl, yl, xl, simdata);
    std::swap(simdata.ntmp, simdata.n);
}

void Streaming::_push_tbb(SimData &simdata) const {
    const auto e = simdata.n->get_extents();
    tbb::parallel_for
    (uint32_t(e.zbegin), e.zend, [this, &e, &simdata] (size_t zl) {
        for (auto yl = e.ybegin; yl < e.yend; ++yl)
            for (auto xl = e.xbegin; xl < e.xend; ++xl)
                _push_kernel(zl, yl, xl, simdata);
    });
    std::swap(simdata.ntmp, simdata.n);
}

void Streaming::_pull_ref(SimData &simdata) const {
    throw std::logic_error("Streaming::_pull_ref has not yet been implemented");
}

void Streaming::_pull_tbb(SimData &simdata) const {
    throw std::logic_error("Streaming::_pull_tbb has not yet been implemented");
}
