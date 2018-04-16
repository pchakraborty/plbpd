#ifndef LBMODEL_HPP
#define LBMODEL_HPP

#include <string>
#include <vector>

class LBModel final{
  
private:
  
    std::string _modelName;
    uint32_t _numVelocityVectors; // 9 for D2Q9, 27 for D3Q27
    float _speedOfSoundSquared;
    std::vector<int32_t> _latticeVelocity; // size=27x3 for D3Q27, 9x2 for D2Q9
    std::vector<float> _directionalWeights; // emacs size=19 for D3Q19
    void setD2Q9();
    // void setD3Q19();
    void setD3Q27();
    // Disable assignment and copy constructors
    
public:
    
    LBModel() = delete;
    LBModel(std::string modelName);
    virtual ~LBModel();
    std::string getModelName() const;
    uint32_t getNumVelocityVectors() const;
    float getSpeedOfSoundSquared() const;
    const std::vector<int32_t> &getLatticeVelocities() const;
    const std::vector<float> &getDirectionalWeights() const;
    
};

#endif
