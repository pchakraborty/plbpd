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

template <typename T, uint32_t num_buffer_layers>
class Field{

private:

    uint32_t _zlen, _ylen, _xlen; // field dimensions
    uint32_t _vlen; // vector length at each field location
    std::vector<T> _arrdata;
    uint32_t _zfactor, _yfactor, _xfactor; // factors for index calculation

    void _initialize(uint32_t zlen, uint32_t ylen, uint32_t xlen, uint32_t vlen, const T& value){
        assert(zlen > 0);
        assert(ylen > 0);
        assert(xlen > 0);
        _zlen = zlen + 2*num_buffer_layers;
        _ylen = ylen + 2*num_buffer_layers;
        _xlen = xlen + 2*num_buffer_layers;
        _vlen = vlen;
        _arrdata.resize(_zlen*_ylen*_xlen*_vlen, value);
        _zfactor = _vlen*_xlen*_ylen;
        _yfactor = _vlen*_xlen;
        _xfactor = _vlen;
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

    inline std::tuple<size_t, size_t, size_t> get_dimensions() const{
        return std::make_tuple(_zlen, _ylen, _xlen);
    }

    inline size_t get_vector_length() const{
        return _vlen;
    }

    inline FieldExtents get_extents() const{
        FieldExtents e;
        // z extents
        e.zbegin = 0 + num_buffer_layers;
        e.zend = _zlen - num_buffer_layers;
        // y extents
        e.ybegin = 0 + num_buffer_layers;
        e.yend = _ylen - num_buffer_layers;
        // x extents
        e.xbegin = 0 + num_buffer_layers;
        e.xend = _xlen - num_buffer_layers;
        return e;
    }

    inline std::tuple<uint32_t, uint32_t, uint32_t>
    get_neighbor(uint32_t z, uint32_t y, uint32_t x, std::array<int, 3> ck) const{
        return std::make_tuple(z+ck[2],y+ck[1],x+ck[0]);
    }

    /* 3D creation and access */

    Field(uint32_t zlen, uint32_t ylen, uint32_t xlen, const T& value){
        _initialize(zlen, ylen, xlen, 1, value);
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x){
        assert(_vlen == 1);
        return _arrdata[z*_zfactor + y*_yfactor + x];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x) const{
        assert(_vlen == 1);
        return _arrdata[z*_zfactor + y*_yfactor + x];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x) const{
        assert(_vlen == 1);
        return &_arrdata[z*_zfactor + y*_yfactor + x];
    }

    inline T* get(uint32_t z, uint32_t y, uint32_t x){
        assert(_vlen == 1);
        return &_arrdata[z*_zfactor + y*_yfactor + x];
    }

    inline uint32_t get_linear_index(uint32_t z, uint32_t y, uint32_t x){
        assert(_vlen == 1);
        return x+(y+z*_ylen)*_xlen;
    }

    /* 4D creation and access */

    // Initialize a field of vectors (of length vlen)
    Field(uint32_t zlen, uint32_t ylen, uint32_t xlen, uint32_t vlen, const T& value){
        assert(vlen > 1);
        _initialize(zlen, ylen, xlen, vlen, value);
    }

    inline T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        assert(_vlen > 1);
        return _arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v];
    }

    inline const T& at(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const{
        assert(_vlen > 1);
        return _arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v];
    }

    inline const T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v) const{
        assert(_vlen > 1);
        return &_arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v];
    }

    inline T* get(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        assert(_vlen > 1);
        return &_arrdata[z*_zfactor + y*_yfactor + x*_xfactor + v];
    }

    inline uint32_t get_linear_index(uint32_t z, uint32_t y, uint32_t x, uint32_t v){
        assert(_vlen > 1);
        return v+(x+(y+z*_ylen)*_xlen)*_vlen;
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
