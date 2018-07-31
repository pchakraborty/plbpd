#include <iostream>
#include <chrono>
#include "Streaming.hpp"
#include "tbb/tbb.h"

float Streaming::_time = 0.0f;

float Streaming::get_total_time() const{
    return _time;
}

Streaming::Streaming(const LBModel *lbmodel, std::string stream_type):
    _lbmodel(lbmodel){
    _reference = false;
    _set_stream_type(stream_type);
}

Streaming::Streaming(const LBModel *lbmodel, std::string stream_type, bool reference):
    _lbmodel(lbmodel), _reference(reference){
    _set_stream_type(stream_type);
}

Streaming::~Streaming(){}

void Streaming::_set_stream_type(std::string stream_type){
    std::vector<std::string> types = {"push", "pull"};
    if (std::find(types.begin(), types.end(), stream_type) != types.end())
        _stream_type = stream_type;
    else
        throw std::invalid_argument("stream_type ["+stream_type+"] not recognized!");
}    

void Streaming::operator()(SimData &simdata) const{
    auto start = std::chrono::system_clock::now();

    if (_reference)
        _streaming_ref(simdata);
    else
        _streaming_tbb(simdata);
    
    std::chrono::duration<float> elapsed = std::chrono::system_clock::now()-start;
    Streaming::_time += elapsed.count();
}

void Streaming::_streaming_ref(SimData &simdata) const{
    if (_stream_type=="pull")
        _pull_ref(simdata);
    else if (_stream_type=="push")
        _push_ref(simdata);
    else
        throw std::logic_error("Streaming::_streaming_ref: unknown stream type "+_stream_type);
}

void Streaming::_streaming_tbb(SimData &simdata) const{
    if (_stream_type=="pull")
        _pull_tbb(simdata);
    else if (_stream_type=="push")
        _push_tbb(simdata);
    else
        throw std::logic_error("Streaming::_streaming_tbb: unknown stream type "+_stream_type);
}

void Streaming::_push_ref(SimData &simdata) const{
    const auto kdim = simdata.n->get_vector_length();
    const auto c = _lbmodel->get_directional_velocities();
    const auto e = simdata.n->get_extents();
    for (auto zl=e.zbegin; zl<e.zend; ++zl){
        for (auto yl=e.ybegin; yl<e.yend; ++yl){
            for (auto xl=e.xbegin; xl<e.xend; ++xl){
                for (auto k=0; k<kdim; ++k){
                    auto ck = &c[k*3];
                    size_t nz, ny, nx;
                    std::tie(nz, ny, nx) = simdata.n->get_neighbor(zl,yl,xl, ck);
                    simdata.ntmp->at(nz,ny,nx,k) = simdata.n->at(zl,yl,xl,k);
                }
            }
        }
    }
    std::swap(simdata.ntmp, simdata.n);
}

void Streaming::_push_tbb(SimData &simdata) const{
    const auto kdim = simdata.n->get_vector_length();
    const auto c = _lbmodel->get_directional_velocities();
    const auto e = simdata.n->get_extents();
    tbb::parallel_for
        (uint32_t(e.zbegin), e.zend, [this, &e, kdim, &c, &simdata] (size_t zl){
            for (auto yl=e.ybegin; yl<e.yend; ++yl){
                for (auto xl=e.xbegin; xl<e.xend; ++xl){
                    for (auto k=0; k<kdim; ++k){
                        auto ck = &c[k*3];
                        size_t nz, ny, nx;
                        std::tie(nz, ny, nx) = simdata.n->get_neighbor(zl,yl,xl, ck);
                        simdata.ntmp->at(nz,ny,nx,k) = simdata.n->at(zl,yl,xl,k);
                    }
                }
            }
        });
    std::swap(simdata.ntmp, simdata.n);
}

void Streaming::_pull_ref(SimData &simdata) const{
    throw std::logic_error("Streaming::_pull_ref has not yet been implemented");
}

void Streaming::_pull_tbb(SimData &simdata) const{
    throw std::logic_error("Streaming::_pull_tbb has not yet been implemented");
}
