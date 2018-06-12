#ifndef FIELD_HPP
#define FIELD_HPP

#include <tuple>
#include <vector>
#include <cassert>
#include <iostream>

namespace Field{

struct FieldExtents{
    uint32_t zbegin, zend;
    uint32_t ybegin, yend;
    uint32_t xbegin, xend;
};

template <typename T, uint32_t _num_buffer_layers>
class Field{

private:

    uint32_t _zlen, _ylen, _xlen; // field dimensions
    uint32_t _vlen; // vector length at each field location
    //uint32_t _num_buffer_layers;
    std::vector<T> _arrdata;
    uint32_t _zfactor, _yfactor, _xfactor, _vfactor; // factors for index calculation

    void _initialize(uint32_t zlen, uint32_t ylen, uint32_t xlen, uint32_t vlen, const T& value){
        // _num_buffer_layers = n_buffer_layers;
        _zlen = zlen + 2*_num_buffer_layers;
        _ylen = ylen + 2*_num_buffer_layers;
        _xlen = xlen + 2*_num_buffer_layers;
        _vlen = vlen;
        _arrdata.resize(_zlen*_ylen*_xlen*_vlen, value);
        _zfactor = _vlen*_xlen*_ylen;
        _yfactor = _vlen*_xlen;
        _xfactor = _vlen;
        _vfactor = 1;
    }

public:

    ~Field(){}

    Field(Field&) = delete;

    Field& operator=(Field&) = delete;

    void print() const{
        for (auto val: _arrdata)
            std::cout<<val<<" ";
        std::cout<<std::endl;
    }

    std::tuple<size_t, size_t, size_t> get_dimensions() const{
        return std::make_tuple(_zlen, _ylen, _xlen);
    }

    size_t get_vector_length() const{
        return _vlen;
    }

    FieldExtents get_extents() const{
        FieldExtents e;
        // z extents
        e.zbegin = 0 + _num_buffer_layers;
        e.zend = _zlen - _num_buffer_layers;
        // y extents
        e.ybegin = 0 + _num_buffer_layers;
        e.yend = _ylen - _num_buffer_layers;
        // x extents
        e.xbegin = 0 + _num_buffer_layers;
        e.xend = _xlen - _num_buffer_layers;
        return e;
    }

    /* 3D creation and access */

    Field(uint32_t zlen, uint32_t ylen, uint32_t xlen, const T& value){
        _initialize(zlen, ylen, xlen, 1, value);
    }
    
    inline T& at(uint32_t z, uint32_t y, uint32_t x){
        assert(_vlen==1);
        return _arrdata[z*_zfactor + y*_yfactor + x*_xfactor];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x) const{
        assert(_vlen==1);
        return _arrdata[z*_zfactor + y*_yfactor + x*_xfactor];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x) const{
        assert(_vlen==1);
        return &_arrdata[z*_zfactor + y*_yfactor + x*_xfactor];
    }

    T* get(uint32_t z, uint32_t y, uint32_t x){
        assert(_vlen==1);
        return &_arrdata[z*_zfactor + y*_yfactor + x*_xfactor];
    }

    /* 4D creation and access */

    // Initialize a 3D (z,y,x) field vectors (of length vlen)
    // vlen=1 => a 3D field of scalars, e.g. density
    // vlen=n => a 3D field of length n vectors
    Field(uint32_t zlen, uint32_t ylen, uint32_t xlen, uint32_t vlen, const T& value){
        assert(vlen>1);
        _initialize(zlen, ylen, xlen, vlen, value);
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        assert(_vlen>1);
        return _arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v*_vfactor];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const{
        assert(_vlen>1);
        return _arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v*_vfactor];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const{
        assert(_vlen>1);
        return &_arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v*_vfactor];
    }
    
    inline T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        assert(_vlen>1);
        return &_arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v*_vfactor];
    }
    
    // Iterate over the vector in Field class
    class iterator {
    private:
        T *_val;
    public:
        inline iterator(T *val) : _val(val) {}
        inline iterator& operator++(){ // pre-fix
            ++_val;
            return *this;
        }
        inline bool operator!=(const iterator& rhs){
            return _val != rhs._val;
        }
        inline T& operator*(){
            return *_val;
        }
    };

    iterator begin(){
        return iterator(_arrdata.data());
    }

    iterator end(){
        return iterator(_arrdata.data()+_zlen*_ylen*_xlen*_vlen);
    }

}; // class Field

} // namespace Field

#endif
