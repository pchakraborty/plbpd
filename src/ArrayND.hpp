#include <vector>
#include <tuple>
#include <stdexcept>
#include <mm_malloc.h>
#include <immintrin.h>

#ifndef ARRAYND_HPP
#define ARRAYND_HPP

namespace ArrayND{

    // 2D array
    template <class T>
    class Array2D{
    private:
        uint32_t _dim1;
        uint32_t _dim2;
        std::vector<T> _arrdata;
    public:
        inline std::tuple<size_t, size_t> getDimensions() const{
            return std::tie(_dim1, _dim2);
        }
        __attribute__((always_inline))
        T& at(int x1, int x2) {
            return _arrdata[x1*_dim2 + x2];
        }
        __attribute__((always_inline))
        const T& at(int x1, int x2) const{
            return _arrdata[x1*_dim2 + x2];
        }
        inline const T* get() const{
            return _arrdata.data();
        }
        Array2D(uint32_t dim1, uint32_t dim2){
            _dim1=dim1;
            _dim2=dim2;
            _arrdata.resize(_dim2*_dim1, static_cast<T>(0.0f));
        }
        ~Array2D(){}
        Array2D(Array2D&) = delete;
        Array2D& operator=(Array2D&) = delete;
    };

    template <class T>
    class Array3D{
    private:
        uint32_t _dim1;
        uint32_t _dim2;
        uint32_t _dim3;
        std::vector<T> _arrdata;
        uint32_t _o1, _o2, _o3;
    public:
        inline std::tuple<size_t, size_t, size_t> getDimensions() const{
            return std::tie(_dim1, _dim2, _dim3);
        }
        __attribute__((always_inline))
        T& at(uint32_t x1, uint32_t x2, uint32_t x3){
            return _arrdata[x1*_o1 + x2*_o2 + x3*_o3];
        }
        __attribute__((always_inline))
        const T& at(uint32_t x1, uint32_t x2, uint32_t x3) const{
            return _arrdata[x1*_o1 + x2*_o2 + x3*_o3];
        }
        inline const T* get() const{
            return _arrdata.data();
        }
        Array3D(uint32_t  dim1, uint32_t dim2, uint32_t dim3){
            _dim1 = dim1;
            _dim2 = dim2;
            _dim3 = dim3;
            _o1=_dim2*_dim3;
            _o2=_dim3;
            _o3=1;
            _arrdata.resize(_dim1*_dim2*_dim3, static_cast<T>(0.0f));
        }
        ~Array3D(){}
        Array3D(Array3D&) = delete;
        Array3D& operator=(Array3D &) = delete;
    };

    template <class T>
    class Array4D{
    private:
        uint32_t _dim1, _dim2, _dim3, _dim4;
        std::vector<T> _arrdata;
        uint32_t _o1, _o2, _o3, _o4;
    public:
        inline std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> getDimensions() const{
            return std::tie(_dim1, _dim2, _dim3, _dim4);
        }
        __attribute__((always_inline))
        T& at(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4){
            return _arrdata[x1*_o1 + x2*_o2 + x3*_o3 + x4*_o4];
        }
        __attribute__((always_inline))
        const T& at(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4) const{
            return _arrdata[x1*_o1 + x2*_o2 + x3*_o3 + x4*_o4];
        }
        inline const T* get() const{
            return _arrdata.data();
        }
        Array4D(uint32_t dim1, uint32_t dim2, uint32_t dim3, uint32_t dim4){
            _dim1=dim1;
            _dim2=dim2;
            _dim3=dim3;
            _dim4=dim4;
            _arrdata.resize(_dim1*_dim2*_dim3*_dim4, static_cast<T>(0.0f));
            _o1=_dim2*_dim3*_dim4 ;
            _o2=_dim3*_dim4 ;
            _o3=_dim4 ;
            _o4=1;
        }
        ~Array4D(){}
        Array4D(Array4D&) = delete;
        Array4D& operator=(Array4D&) = delete;
    };

    template <class T>
    class Array4D_simd{
    private:
        uint32_t _dim1, _dim2, _dim3, _dim4;
        T* _arrdata;
        __m128i _oi; // 4 32-bit ints
        __attribute__((always_inline))
        int32_t _xioi(const __m128i _xi) const{
            auto _pi = _mm_mullo_epi32(_xi, _oi);
            // auto _ri = _mm_hadd_epi32(_pi, _pi);
            // _ri = _mm_hadd_epi32(_ri, _ri);
            auto _ri = _mm_add_epi32(_pi, _mm_srli_si128(_pi, 8));
            _ri = _mm_add_epi32(_ri, _mm_srli_si128(_ri, 4));
            return _mm_extract_epi32(_ri, 0);
        }
    public:
        inline std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> getDimensions() const{
            return std::tie(_dim1, _dim2, _dim3, _dim4);
        }
       __attribute__((always_inline))
       T& at(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4){
           auto xi = _mm_set_epi32(x4, x3, x2, x1);
           return _arrdata[_xioi(_mm_set_epi32(x4, x3, x2, x1))];
       }
        __attribute__((always_inline))
        const T& at(uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4) const{
            auto xi = _mm_set_epi32(x4, x3, x2, x1);
            auto ndx = _xioi(xi);
            return _arrdata[ndx];
        }
        inline const T* get() const{
            return _arrdata;
        }
        Array4D_simd(uint32_t dim1, uint32_t dim2, uint32_t dim3, uint32_t dim4){
            _dim1=dim1;
            _dim2=dim2;
            _dim3=dim3;
            _dim4=dim4;
            const auto alignment = 64;
            auto arrsize = _dim1*_dim2*_dim3*_dim4*sizeof(T);
            if (arrsize%alignment!=0)
                throw std::invalid_argument("Array4D cannot be aligned to 64 byte bdry");
            _arrdata = static_cast<T*>(_mm_malloc(arrsize, alignment));
            _oi = _mm_set_epi32(1, _dim4, _dim3*_dim4, _dim2*_dim3*_dim4);
        }
        ~Array4D_simd(){
            if (_arrdata) _mm_free(_arrdata);
        }
        Array4D_simd(Array4D_simd&) = delete;
        Array4D_simd& operator=(Array4D_simd&) = delete;
    };

} // namespace ArrayND

#endif

