#ifndef FIELD_HPP
#define FIELD_HPP

#include <tuple>
#include <vector>
#include <cassert>
#include <iostream>
#include <stdexcept>

namespace Field{

struct FieldExtents{
    uint32_t zbegin, zend;
    uint32_t ybegin, yend;
    uint32_t xbegin, xend;
};

template <typename T, uint32_t num_buffer_layers>
class ScalarField{

private:

    uint32_t _zlen, _ylen, _xlen; // field dimensions
    std::vector<T> _arrdata;
    uint32_t _zfactor, _yfactor; // factors for index calculation
    FieldExtents _e;

public:

    ScalarField(uint32_t zlen, uint32_t ylen, uint32_t xlen, const T& value){

        if (xlen == 0 || ylen==0 || zlen ==0){
            throw std::invalid_argument("at least one of xlen/ylen/zlen is zero");
        }

        _zlen = zlen + 2*num_buffer_layers;
        _ylen = ylen + 2*num_buffer_layers;
        _xlen = xlen + 2*num_buffer_layers;

        _arrdata.resize(_zlen*_ylen*_xlen, value);

        _zfactor = _xlen*_ylen;
        _yfactor = _xlen;

        _e.zbegin = 0 + num_buffer_layers;
        _e.zend = _zlen - num_buffer_layers;
        _e.ybegin = 0 + num_buffer_layers;
        _e.yend = _ylen - num_buffer_layers;
        _e.xbegin = 0 + num_buffer_layers;
        _e.xend = _xlen - num_buffer_layers;
    }

    ~ScalarField(){}

    ScalarField(ScalarField&) = delete;

    ScalarField& operator=(ScalarField&) = delete;

    inline std::tuple<uint32_t, uint32_t, uint32_t> get_dimensions() const{
        return std::make_tuple(_zlen, _ylen, _xlen);
    }

    inline const FieldExtents& get_extents() const{
        return _e;
    }

    inline uint32_t sub2ind(uint32_t z, uint32_t y, uint32_t x){
        return z*_zfactor + y*_yfactor + x;
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x){
        return _arrdata[sub2ind(z,y,x)];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x) const{
        return _arrdata[sub2ind(z,y,x)];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x) const{
        return &_arrdata[sub2ind(z,y,x)];
    }

    inline T* get(uint32_t z, uint32_t y, uint32_t x){
        return &_arrdata[sub2ind(z,y,x)];
    }

}; // class SclalarField

template <typename T, uint32_t num_buffer_layers>
class VectorField{

private:

    uint32_t _zlen, _ylen, _xlen; // field dimensions
    uint32_t _vlen; // vector length at each field location
    std::vector<T> _arrdata;
    uint32_t _zfactor, _yfactor, _xfactor; // factors for index calculation
    FieldExtents _e;
    
    void _check_arg(uint32_t arg){
        if (arg == 0){
            throw std::invalid_argument("argument is zero");
        }
    }
    
    void _set_lengths(uint32_t zlen, uint32_t ylen, uint32_t xlen, uint32_t vlen){
        _zlen = zlen + 2*num_buffer_layers;
        _ylen = ylen + 2*num_buffer_layers;
        _xlen = xlen + 2*num_buffer_layers;
        _vlen = vlen;
    }
    
    void _set_factors(){
        _zfactor = _vlen*_xlen*_ylen;
        _yfactor = _vlen*_xlen;
        _xfactor = _vlen;
    }        

    void _set_extents(){
        _e.zbegin = 0 + num_buffer_layers;
        _e.zend = _zlen - num_buffer_layers;
        _e.ybegin = 0 + num_buffer_layers;
        _e.yend = _ylen - num_buffer_layers;
        _e.xbegin = 0 + num_buffer_layers;
        _e.xend = _xlen - num_buffer_layers;
    }        
    
public:

    // Initialize a field of vectors (of length vlen)
    VectorField(uint32_t zlen, uint32_t ylen, uint32_t xlen, uint32_t vlen, const T& value){
        for (auto arg: std::vector<uint32_t>{zlen, ylen, xlen, vlen})
            _check_arg(arg);
        _set_lengths(zlen, ylen, xlen, vlen);
        _arrdata.resize(_zlen*_ylen*_xlen*_vlen, value);
        _set_factors();
        _set_extents();
    }

    ~VectorField(){}

    VectorField(VectorField&) = delete;

    VectorField& operator=(VectorField&) = delete;

    inline std::tuple<uint32_t, uint32_t, uint32_t> get_dimensions() const{
        return std::make_tuple(_zlen, _ylen, _xlen);
    }

    inline uint32_t get_vector_length() const{
        return _vlen;
    }

    inline const FieldExtents& get_extents() const{
        return _e;
    }

    inline uint32_t sub2ind(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        return z*_zfactor + y*_yfactor + x*_xfactor + v;
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        return _arrdata[sub2ind(z,y,x,v)];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const{
        return _arrdata[sub2ind(z,y,x,v)];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const{
        return &_arrdata[sub2ind(z,y,x,v)];
    }

    inline T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        return &_arrdata[sub2ind(z,y,x,v)];
    }

    inline std::tuple<uint32_t, uint32_t, uint32_t>
    get_neighbor(uint32_t z, uint32_t y, uint32_t x, const int* ck) const{
        return std::make_tuple(z+ck[2],y+ck[1],x+ck[0]);
    }

}; // class VectorField

} // namespace Field

#endif
