#include "LBModel.hpp"
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
    _numVelocityVectors = 9;
    _speedOfSoundSquared = 1./3.;
    _latticeVelocity = {
         0.,  0., 0.,  //  0
         1.,  0., 0.,  //  1
        -1.,  0., 0.,  //  2
         0.,  1., 0.,  //  3
         0., -1., 0.,  //  4
         1.,  1., 0.,  //  5
        -1., -1., 0.,  //  6
         1., -1., 0.,  //  7
        -1.,  1., 0.   //  8
    };
    _directionalWeights = {
        4./9.,
        1./9.,  1./9.,  1./9.,  1./9.,
        1./36., 1./36., 1./36., 1./36.
    };
}

void LBModel::setD3Q27(){
    _numVelocityVectors = 27;
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
    _directionalWeights = {
        8./27., // rest particle
        // group i
        2./27.,  2./27.,  2./27.,  2./27., 2./27., 2./27.,
        // group ii
        1./54.,  1./54.,  1./54.,  1./54., 1./54., 1./54.,
        1./54.,  1./54.,  1./54.,  1./54., 1./54., 1./54.,
        // group iii
        1./216., 1./216., 1./216., 1./216.,
        1./216., 1./216., 1./216., 1./216.
    };
}

std::string LBModel::getModelName() const{
    return _modelName;
}

unsigned short LBModel::getNumVelocityVectors() const{
    return _numVelocityVectors;
}

double LBModel::getSpeedOfSoundSquared() const{
    return _speedOfSoundSquared;
}

const std::vector<double> &LBModel::getLatticeVelocities() const{
    return _latticeVelocity;
}

const std::vector<double> &LBModel::getDirectionalWeights() const{
    return _directionalWeights;
}
