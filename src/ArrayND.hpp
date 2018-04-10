#include <cassert>
#include <vector>
#include <tuple>

#ifndef ARRAYND_HPP
#define ARRAYND_HPP

namespace ArrayND{

    // 2D array
    template <class T>
    class Array2D{
    private:
        int _dim1;
        int _dim2;
        std::vector<T> arrdata;
    public:
        inline std::tuple<size_t, size_t> getDimensions() const{
            return std::tie(_dim1, _dim2);
        }
        __attribute__((always_inline))
        T& at(int x1, int x2) {
            return arrdata[x1*_dim2 + x2];
        }
        __attribute__((always_inline))
        const T& at(int x1, int x2) const{
            return arrdata[x1*_dim2 + x2];
        }
        inline const T* get() const{
            return arrdata.data();
        }
        Array2D(int dim1, int dim2){
            assert(dim1>0);
            assert(dim2>0);
            _dim1=dim1;
            _dim2=dim2;
            arrdata.resize(_dim2*_dim1, static_cast<T>(0.0f));
        }
        ~Array2D(){}
        Array2D(Array2D&) = delete;
        Array2D& operator=(Array2D&) = delete;
    };

    template <class T>
    class Array3D{
    private:
        int _dim1;
        int _dim2;
        int _dim3;
        std::vector<T> arrdata;
        size_t _o1, _o2, _o3;
    public:
        inline std::tuple<size_t, size_t, size_t> getDimensions() const{
            return std::tie(_dim1, _dim2, _dim3);
        }
        __attribute__((always_inline))
        T& at(int x1, int x2, int x3){
            return arrdata[x1*_o1 + x2*_o2 + x3*_o3];
        }
        __attribute__((always_inline))
        const T& at(int x1, int x2, int x3) const{
            return arrdata[x1*_o1 + x2*_o2 + x3*_o3];
        }
        inline const T* get() const{
            return arrdata.data();
        }
        Array3D(int  dim1, int dim2, int dim3){
            assert(dim1>0);
            assert(dim2>0);
            assert(dim3>0);
            _dim1 = dim1;
            _dim2 = dim2;
            _dim3 = dim3;
            _o1=_dim2*_dim3;
            _o2=_dim3;
            _o3=1;
            arrdata.resize(_dim1*_dim2*_dim3, static_cast<T>(0.0f));
        }
        ~Array3D(){}
        Array3D(Array3D&) = delete;
        Array3D& operator=(Array3D &) = delete;
    };

    template <class T>
    class Array4D{
    private:
        int _dim1;
        int _dim2;
        int _dim3;
        int _dim4;
        std::vector<T> arrdata;
        size_t _o1, _o2, _o3, _o4;
    public:
        inline std::tuple<size_t, size_t, size_t, size_t> getDimensions() const{
            return std::tie(_dim1, _dim2, _dim3, _dim4);
        }
        __attribute__((always_inline)) T& at(int x1, int x2, int x3, int x4){
            return arrdata[x1*_o1 + x2*_o2 + x3*_o3 + x4*_o4 ];
        }
        __attribute__((always_inline)) const T& at(int x1, int x2, int x3, int x4) const {
            return arrdata[x1*_o1 + x2*_o2 + x3*_o3 + x4*_o4 ];
        }
        inline const T* get() const{
            return arrdata.data();
        }
        Array4D(int  dim1, int dim2, int dim3, int dim4){
            assert(dim1>0);
            assert(dim2>0);
            assert(dim3>0);
            assert(dim4>0);
            _dim1=dim1;
            _dim2=dim2;
            _dim3=dim3;
            _dim4=dim4;
            arrdata.resize(_dim1*_dim2*_dim3*_dim4, static_cast<T>(0.0f));
            _o1=_dim2*_dim3*_dim4 ;
            _o2=_dim3*_dim4 ;
            _o3=_dim4 ;
            _o4=1;
        }
        ~Array4D(){}
        Array4D(Array4D&) = delete;
        Array4D& operator=(Array4D&) = delete;
    };

} // namespace ArrayND

#endif

