#ifndef SRC_BLOCK3D_HPP_
#define SRC_BLOCK3D_HPP_

#include <tuple>
#include <string>

class Block3d final {

private:
    uint32_t _xstart, _ystart, _zstart; // block origin
    uint32_t _xsize, _ysize, _zsize;  // block size in each direction

public:
    Block3d() = delete;
    Block3d(size_t xdim, size_t ydim, size_t zdim);
    Block3d(Block3d&) = delete;
    Block3d& operator=(Block3d&) = delete;
    ~Block3d();
    std::tuple<size_t, size_t, size_t> get_dimensions() const;
};

#endif  // SRC_BLOCK3D_HPP_
