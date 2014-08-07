#ifndef TIMER_HPP
#define TIMER_HPP

#include "pointers.hpp"

// TIME_N should always be the last term. it defines the array size
enum{MAINLOOP, COMM, COLLIDENSTREAM_FL, APPLYBC_FL, CALCMOMENTS_FL,
     COLLIDENSTREAM_TH, APPLYBC_TH, CALCMOMENTS_TH, OUTPUT, 
     INTEGRATE, RESET_FLAG_MOMENT, TIME_N};

class Timer: protected Pointers{
public:
  double *elapsed;

  Timer(PLBPD* plbpd);
  ~Timer();

  void setup();
  void start(int);
  void end(int);

private:
  double *_start_, _end_;
};

#endif
