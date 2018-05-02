#ifndef LBMODEL_HPP
#define LBMODEL_HPP

#include <string>
#include <vector>

class LBModel final{
  
private:
  
    std::string _modelName;
    uint32_t _numberOfDirections; // 9 for D2Q9, 27 for D3Q27
    uint32_t _numberOfDirectionsSse; // 28 for D3Q27
    uint32_t _numberOfDirectionsAvx; // 32 for D3Q27
    float _speedOfSoundSquared;
    std::vector<int32_t> _latticeVelocity; // size=27x3 for D3Q27, 9x2 for D2Q9
    std::vector<float> _latticeVelocitySse; // size=27*4 for D3Q27 with padding for SSE
    std::vector<float> _directionalWeights; // size=19 for D3Q19
    std::vector<float> _directionalWeightsAvx; // size=32 for D3Q27 (zero-padded for AVX)
    std::vector<uint32_t> _reverse;
    void setD2Q9();
    void setD3Q19();
    void setD3Q27();
    // Disable assignment and copy constructors
    
public:
    
    LBModel() = delete;
    LBModel(std::string modelName);
    virtual ~LBModel();
    std::string getModelName() const;
    uint32_t getNumberOfDirections() const;
    uint32_t getNumberOfDirectionsSse() const;
    uint32_t getNumberOfDirectionsAvx() const;
    float getSpeedOfSoundSquared() const;
    const std::vector<int32_t> &getLatticeVelocities() const;
    const std::vector<float> &getLatticeVelocitiesSse() const;
    const std::vector<float> &getDirectionalWeights() const;
    const std::vector<float> &getDirectionalWeightsAvx() const;
    const std::vector<uint32_t> &getReverse() const;
};

#endif
