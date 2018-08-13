#include <cassert>
#include "tbb/tbb.h"
#include "CalcMoments.hpp"

float CalcMoments::_time_calc_moment = 0.0f;

CalcMoments::CalcMoments(const LBModel *lbmodel): _lbmodel(lbmodel) {}

CalcMoments::~CalcMoments() {}

float CalcMoments::get_total_time() const {
    return _time_calc_moment;
}
