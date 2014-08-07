#ifndef COMM_HPP
#define COMM_HPP

#include "mpi.h"
#include "pointers.hpp"

class Comm: protected Pointers{
public:
  MPI_Comm world;                // communicator for all processes
  MPI_Info info;
  int me, nprocs;                // my place among all procs
  bool root;                     // whether I am root or not
  int npx, npy, npz;             // user supplied processor grid
  int myloc[3];                  // which proc I am in, in each dim
  int east, west;                // my 6 nbring procs
  int north, south;
  int up, down;

  Comm(PLBPD *, MPI_Comm, MPI_Info);
  ~Comm();

  void setup();                  // setup 3D grid of procs
  void alloc();                  // allocate memory for send and recv buffers
  void communicate();            // exchange halo region between procs
  void communicate_serial();

private:
  int maxbufsize;                // max size of send and recv buffers
  int nDoubles;                  // num dbles at each node
  // send and recv buffers
  int count_f, count_b;
  double *sendbuf_f, *sendbuf_b, *recvbuf_f, *recvbuf_b;
  void pack(const int direction);
  void unpack(const int direction);
};

#endif
