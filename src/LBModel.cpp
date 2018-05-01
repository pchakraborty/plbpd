#include "LBModel.hpp"
#include <mm_malloc.h>
#include <stdexcept>

LBModel::LBModel(std::string myModelName){
    _modelName = myModelName;
    if (_modelName=="D2Q9"){
        setD2Q9();
    } else if (_modelName=="D3Q27"){
        setD3Q27();
    } else{
        throw std::logic_error("Model "+_modelName+" has not been implemented");
    }
}

LBModel::~LBModel(){}

void LBModel::setD2Q9(){
    _numberOfDirections = 9;
    _speedOfSoundSquared = 1./3.;
    _latticeVelocity = { // y-component is 0
         0, 0, 0,  //  0
         1, 0, 0,  //  1
        -1, 0, 0,  //  2
         0, 0, 1,  //  3
         0, 0,-1,  //  4
         1, 0, 1,  //  5
        -1, 0,-1,  //  6
         1, 0,-1,  //  7
        -1, 0, 1   //  8
    };
    _directionalWeights = {
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

void LBModel::setD3Q27(){
    _numberOfDirections = 27;
    _numberOfDirectionsSse = 28;
    _numberOfDirectionsAvx = 32;
    _speedOfSoundSquared = 1./3.;

    _latticeVelocity = {
         0, 0, 0,   //  0
        // group i
         1, 0, 0,   //  1
        -1, 0, 0,   //  2
         0, 1, 0,   //  3
         0,-1, 0,   //  4
         0, 0, 1,   //  5
         0, 0,-1,   //  6
        // group ii
         1, 1, 0,   //  7
        -1,-1, 0,   //  8
         1,-1, 0,   //  9
        -1, 1, 0,   // 10
         1, 0, 1,   // 11
        -1, 0,-1,   // 12
         1, 0,-1,   // 13
        -1, 0, 1,   // 14
         0, 1, 1,   // 15
         0,-1,-1,   // 16
         0, 1,-1,   // 17
         0,-1, 1,   // 18
        // group iii
         1, 1, 1,   // 19
        -1,-1,-1,   // 20
         1, 1,-1,   // 21
        -1,-1, 1,   // 22
         1,-1, 1,   // 23
        -1, 1,-1,   // 24
        -1, 1, 1,   // 25
         1,-1,-1    // 26
    };

    _latticeVelocitySse = {
         0, 0, 0, 0, //  0
        // group i
         1, 0, 0, 0, //  1
        -1, 0, 0, 0, //  2
         0, 1, 0, 0, //  3
         0,-1, 0, 0, //  4
         0, 0, 1, 0, //  5
         0, 0,-1, 0, //  6
        // group ii
         1, 1, 0, 0, //  7
        -1,-1, 0, 0, //  8
         1,-1, 0, 0, //  9
        -1, 1, 0, 0, // 10
         1, 0, 1, 0, // 11
        -1, 0,-1, 0, // 12
         1, 0,-1, 0, // 13
        -1, 0, 1, 0, // 14
         0, 1, 1, 0, // 15
         0,-1,-1, 0, // 16
         0, 1,-1, 0, // 17
         0,-1, 1, 0, // 18
        // group iii
         1, 1, 1, 0, // 19
        -1,-1,-1, 0, // 20
         1, 1,-1, 0, // 21
        -1,-1, 1, 0, // 22
         1,-1, 1, 0, // 23
        -1, 1,-1, 0, // 24
        -1, 1, 1, 0, // 25
         1,-1,-1, 0, // 26
         0, 0, 0, 0  // 27
    };
    
    _directionalWeights = {
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

    _directionalWeightsAvx = {
        // rest particle
        8./27.,
        // group i
        2./27.,  2./27.,  2./27.,  2./27., 2./27., 2./27.,
        // group ii
        1./54.,  1./54.,  1./54.,  1./54., 1./54., 1./54.,
        1./54.,  1./54.,  1./54.,  1./54., 1./54., 1./54.,
        // group iii
        1./216., 1./216., 1./216., 1./216.,
        1./216., 1./216., 1./216., 1./216.,
        // zero padding for avx
        0.0, 0.0, 0.0, 0.0, 0.0
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

std::string LBModel::getModelName() const{
    return _modelName;
}

uint32_t LBModel::getNumberOfDirections() const{
    return _numberOfDirections;
}

uint32_t LBModel::getNumberOfDirectionsSse() const{
    return _numberOfDirectionsSse;
}

uint32_t LBModel::getNumberOfDirectionsAvx() const{
    return _numberOfDirectionsAvx;
}

float LBModel::getSpeedOfSoundSquared() const{
    return _speedOfSoundSquared;
}

const std::vector<int32_t> &LBModel::getLatticeVelocities() const{
    return _latticeVelocity;
}

const std::vector<float> &LBModel::getLatticeVelocitiesSse() const{
    return _latticeVelocitySse;
}

const std::vector<float> &LBModel::getDirectionalWeights() const{
    return _directionalWeights;
}

const std::vector<float> &LBModel::getDirectionalWeightsAvx() const{
    return _directionalWeightsAvx;
}

const std::vector<uint32_t> &LBModel::getReverse() const{
    return _reverse;
}
