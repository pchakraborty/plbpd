#ifndef STREAMING_HPP
#define STREAMING_HPP

#include <string>
#include "SimData.hpp"
#include "LBModel.hpp"

class Streaming{
    
private:

    const LBModel *_lbmodel;
    static float _time;
    bool _reference;
    std::string _stream_type; // push/pull
    
    void _set_stream_type(std::string stream_type);
    void _push_ref(SimData &simdata) const;
    void _push_tbb(SimData &simdata) const;
    void _pull_ref(SimData &simdata) const;
    void _pull_tbb(SimData &simdata) const;
    void _streaming_ref(SimData &simdata) const;
    void _streaming_tbb(SimData &simdata) const;

public:
    
    Streaming(const LBModel *lbmodel, std::string stream_type);
    Streaming(const LBModel *lbmodel, std::string stream_type, bool reference);
    Streaming(Streaming&) = delete;
    Streaming& operator=(Streaming&) = delete;
    ~Streaming();
    void operator()(SimData &simdata) const;
    float get_total_time() const;

};

#endif
