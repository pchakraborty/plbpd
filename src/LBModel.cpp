#include "LBModel.hpp"

#include <mm_malloc.h>
#include <stdexcept>
#include <string>
#include <vector>

LBModel::LBModel(std::string model_name) {
    _model_name = model_name;
    if (_model_name == "D2Q9") {
        set_d2q9();
    } else if (_model_name == "D3Q27") {
        set_d3q27();
    } else if (_model_name == "D3Q19") {
        set_d3q19();
    } else {
        throw std::logic_error("Model "+_model_name+" not yet implemented");
    }
}

LBModel::~LBModel() {}

void LBModel::set_d2q9() {
    _num_directions = 9;
    _speed_of_sound_squared = 1./3.;
    _directional_velocities = {  // y-component is 0
        { 0,  0,  0},  //  0
        { 1,  0,  0},  //  1
        {-1,  0,  0},  //  2
        { 0,  0,  1},  //  3
        { 0,  0, -1},  //  4
        { 1,  0,  1},  //  5
        {-1,  0, -1},  //  6
        { 1,  0, -1},  //  7
        {-1,  0,  1}   //  8
    };
    _directional_weights = {
        4./9.,
        1./9.,  1./9.,  1./9.,  1./9.,
        1./36., 1./36., 1./36., 1./36.
    };
    _reverse = {
        0,
        2,  1,
        4,  3,
        6,  5,
        8,  7
    };
}

void LBModel::set_d3q27() {
    _num_directions = 27;
    _speed_of_sound_squared = 1./3.;

    _directional_velocities = {
        { 0,  0,  0},   //  0
        // group i
        { 1,  0,  0},   //  1
        {-1,  0,  0},   //  2
        { 0,  1,  0},   //  3
        { 0, -1,  0},   //  4
        { 0,  0,  1},   //  5
        { 0,  0, -1},   //  6
        // group ii
        { 1,  1,  0},   //  7
        {-1, -1,  0},   //  8
        { 1, -1,  0},   //  9
        {-1,  1,  0},   // 10
        { 1,  0,  1},   // 11
        {-1,  0, -1},   // 12
        { 1,  0, -1},   // 13
        {-1,  0,  1},   // 14
        { 0,  1,  1},   // 15
        { 0, -1, -1},   // 16
        { 0,  1, -1},   // 17
        { 0, -1,  1},   // 18
        // group iii
        { 1,  1,  1},   // 19
        {-1, -1, -1},   // 20
        { 1,  1, -1},   // 21
        {-1, -1,  1},   // 22
        { 1, -1,  1},   // 23
        {-1,  1, -1},   // 24
        {-1,  1,  1},   // 25
        { 1, -1, -1}    // 26
    };

    _directional_weights = {
        // rest particle
        8./27.,
        // group i
        2./27.,  2./27.,  2./27.,  2./27., 2./27., 2./27.,
        // group ii
        1./54.,  1./54.,  1./54.,  1./54., 1./54., 1./54.,
        1./54.,  1./54.,  1./54.,  1./54., 1./54., 1./54.,
        // group iii
        1./216., 1./216., 1./216., 1./216.,
        1./216., 1./216., 1./216., 1./216.
    };

    _reverse = {
        0,
        2,  1,
        4,  3,
        6,  5,
        8,  7,
        10,  9,
        12, 11,
        14, 13,
        16, 15,
        18, 17,
        20, 19,
        22, 21,
        24, 23,
        26, 25
    };
}

void LBModel::set_d3q19() {
    _num_directions = 19;
    _speed_of_sound_squared = 1./3.;

    _directional_velocities = {
        { 0,  0,  0},   //  0
        { 1,  0,  0},   //  1
        {-1,  0,  0},   //  2
        { 0,  1,  0},   //  3
        { 0, -1,  0},   //  4
        { 0,  0,  1},   //  5
        { 0,  0, -1},   //  6
        { 1,  1,  0},   //  7
        {-1, -1,  0},   //  8
        { 1, -1,  0},   //  9
        {-1,  1,  0},   // 10
        { 1,  0,  1},   // 11
        {-1,  0, -1},   // 12
        { 1,  0, -1},   // 13
        {-1,  0,  1},   // 14
        { 0,  1,  1},   // 15
        { 0, -1, -1},   // 16
        { 0,  1, -1},   // 17
        { 0, -1,  1}    // 18
    };

    _directional_weights = {
        1./3.,
        1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
        1./36., 1./36., 1./36., 1./36.,
        1./36., 1./36., 1./36., 1./36.,
        1./36., 1./36., 1./36., 1./36.,
    };

    _reverse = {
        0,
        2,  1,
        4,  3,
        6,  5,
        8,  7,
        10,  9,
        12, 11,
        14, 13,
        16, 15,
        18, 17
    };
}

std::string LBModel::get_model_name() const {
    return _model_name;
}

uint32_t LBModel::get_num_directions() const {
    return _num_directions;
}

float LBModel::get_speed_of_sound_squared() const {
    return _speed_of_sound_squared;
}

const std::vector<std::array<int32_t, 3> > &LBModel::get_directional_velocities() const {
    return _directional_velocities;
}

const std::vector<float> &LBModel::get_directional_weights() const {
    return _directional_weights;
}

const std::vector<uint32_t> &LBModel::get_reverse() const {
    return _reverse;
}
