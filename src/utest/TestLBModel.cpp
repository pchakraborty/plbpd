#include "../LBModel.hpp"
#include <memory>
#include <cassert>
#include <vector>

int main(){
    { // D2Q9
        auto lbmodel = std::make_unique<LBModel>("D2Q9");
        assert(lbmodel->getModelName()=="D2Q9");
        assert(lbmodel->getNumberOfDirections()==9);
        assert(lbmodel->getSpeedOfSoundSquared()==1.0f/3.0f);
        std::vector<int32_t> latticeVelocity = {
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
        std::vector<float> directionalWeights = {
            4./9.,
            1./9.,  1./9.,  1./9.,  1./9.,
            1./36., 1./36., 1./36., 1./36.
        };
        std::vector<uint32_t> reverse = {
            0,
            2,  1,
            4,  3,
            6,  5,
            8,  7
        };
        assert(lbmodel->getLatticeVelocities()==latticeVelocity);
        assert(lbmodel->getDirectionalWeights()==directionalWeights);
        assert(lbmodel->getReverse()==reverse);
    }

    { // D3Q19
        auto lbmodel = std::make_unique<LBModel>("D3Q19");
        assert(lbmodel->getModelName()=="D3Q19");
        assert(lbmodel->getNumberOfDirections()==19);
        assert(lbmodel->getSpeedOfSoundSquared()==1.0f/3.0f);
        std::vector<int32_t> latticeVelocity = {
             0, 0, 0,   //  0
             1, 0, 0,   //  1
            -1, 0, 0,   //  2
             0, 1, 0,   //  3
             0,-1, 0,   //  4
             0, 0, 1,   //  5
             0, 0,-1,   //  6
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
             0,-1, 1  // 18
        };
        std::vector<float> directionalWeights = {
            1./3.,
            1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
            1./36., 1./36., 1./36., 1./36.,
            1./36., 1./36., 1./36., 1./36.,
            1./36., 1./36., 1./36., 1./36.,
        };
        std::vector<uint32_t> reverse = {
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
        assert(lbmodel->getLatticeVelocities()==latticeVelocity);
        assert(lbmodel->getDirectionalWeights()==directionalWeights);
        assert(lbmodel->getReverse()==reverse);
    }
    
    { // D3Q27
        auto lbmodel = std::make_unique<LBModel>("D3Q27");
        assert(lbmodel->getModelName()=="D3Q27");
        assert(lbmodel->getNumberOfDirections()==27);
        assert(lbmodel->getSpeedOfSoundSquared()==1.0f/3.0f);
        std::vector<int32_t> latticeVelocity = {
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
        std::vector<float> directionalWeights = {
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
        std::vector<uint32_t> reverse = {
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
        assert(lbmodel->getLatticeVelocities()==latticeVelocity);
        assert(lbmodel->getDirectionalWeights()==directionalWeights);
        assert(lbmodel->getReverse()==reverse);
    }

    

    return 0;
}
