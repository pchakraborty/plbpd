#include "timer.hpp"

#include "mpi.h"

/******************************************************************************/

Timer::Timer(PLBPD* plbpd): Pointers(plbpd){
  _start_ = NULL;
  elapsed = NULL; 
}

/******************************************************************************/

Timer::~Timer(){
  if(_start_) delete [] _start_;
  if(elapsed) delete [] elapsed;
}

/******************************************************************************/

void Timer::setup(){
  _start_ = new double[TIME_N];
  elapsed = new double[TIME_N];
  for(int i=0; i<TIME_N; i++){
    _start_[i] = 0.0;
    elapsed[i] = 0.0;
  }
}

/******************************************************************************/

void Timer::start(int which){
  _start_[which] = MPI_Wtime();
}

/******************************************************************************/

void Timer::end(int which){
  _end_ = MPI_Wtime();
  elapsed[which] += (_end_ - _start_[which]);
}

/******************************************************************************/
