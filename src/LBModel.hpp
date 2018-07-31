#ifndef LBMODEL_HPP
#define LBMODEL_HPP

#include <string>
#include <vector>

class LBModel final{
  
private:
  
    std::string _model_name;
    uint32_t _num_directions; // 9 for D2Q9, 27 for D3Q27
    float _speed_of_sound_squared;
    // A nicer representation would be
    //   std::vector<std::array<int32_t, 3> > _lattice_velocity
    // but that is slower than the current implementation
    std::vector<int32_t> _directional_velocities; // size=27x3 for D3Q27, 9x2 for D2Q9
    std::vector<float> _directional_weights; // size=19 for D3Q19
    std::vector<uint32_t> _reverse;
    void set_d2q9();
    void set_d3q19();
    void set_d3q27();
    // Disable assignment and copy constructors
    
public:
    
    LBModel() = delete;
    LBModel(std::string model_name);
    LBModel(LBModel&) = delete;
    LBModel& operator=(LBModel&) = delete;
    virtual ~LBModel();
    std::string get_model_name() const;
    uint32_t get_num_directions() const;
    float get_speed_of_sound_squared() const;
    const std::vector<int32_t> &get_directional_velocities() const;
    const std::vector<float> &get_directional_weights() const;
    const std::vector<uint32_t> &get_reverse() const;
};

#endif
