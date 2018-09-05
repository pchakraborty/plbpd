#ifndef SRC_LBMODEL_HPP_
#define SRC_LBMODEL_HPP_

#include <string>
#include <vector>
#include <array>

class LBModel final {
 private:
    std::string _model_name;
    uint32_t _num_directions;  // 9 for D2Q9, 27 for D3Q27
    float _speed_of_sound_squared;
    // A faster representation would be
    // std::vector<int32_t> _lattice_velocity; where size=19*3 for D3Q19
    // but the current representation makes for more readable code
    std::vector<std::array<int32_t, 3> > _directional_velocities;
    std::vector<float> _directional_weights;  // size=19 for D3Q19
    std::vector<uint32_t> _reverse;
    void set_d2q9();
    void set_d3q19();
    void set_d3q27();

 public:
    LBModel() = delete;
    explicit LBModel(std::string model_name);
    LBModel(LBModel&) = delete;
    LBModel& operator=(LBModel&) = delete;
    virtual ~LBModel();
    std::string get_model_name() const;
    uint32_t get_num_directions() const;
    float get_speed_of_sound_squared() const;
    const std::vector<std::array<int32_t, 3> > &get_directional_velocities() const;
    const std::vector<float> &get_directional_weights() const;
    const std::vector<uint32_t> &get_reverse() const;
};

#endif  // SRC_LBMODEL_HPP_
