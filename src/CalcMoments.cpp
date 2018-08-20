#include <cassert>
#include "tbb/tbb.h"
#include "CalcMoments.hpp"

float CalcMoments::_time_calc_moment = 0.0f;

CalcMoments::CalcMoments(const LBModel *lbmodel)
    : _c (lbmodel->get_directional_velocities()),
      _w (lbmodel->get_directional_weights()),
      _kdim (lbmodel->get_num_directions()) {}

CalcMoments::~CalcMoments() {}

float CalcMoments::get_total_time() const {
    return _time_calc_moment;
}
