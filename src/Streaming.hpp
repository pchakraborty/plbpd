#ifndef STREAMING_HPP
#define STREAMING_HPP

#include <string>
#include "Lattice.hpp"

class Streaming{
    
private:

    const LBModel *_lbmodel;
    static float _time;
    bool _reference;
    std::string _stream_type; // push/pull
    
    void _set_stream_type(std::string stream_type);
    void _push_ref(Lattice &lattice) const;
    void _push_tbb(Lattice &lattice) const;
    void _pull_ref(Lattice &lattice) const;
    void _pull_tbb(Lattice &lattice) const;
    void _streaming_ref(Lattice &lattice) const;
    void _streaming_tbb(Lattice &lattice) const;

public:
    
    Streaming(const LBModel *lbmodel, std::string stream_type);
    Streaming(const LBModel *lbmodel, std::string stream_type, bool reference);
    Streaming& operator=(Streaming&) = delete;
    ~Streaming();
    void operator()(Lattice &lattice) const;
    float get_total_time() const;

};

#endif
