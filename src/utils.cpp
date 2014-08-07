#include "utils.hpp"
#include <cstdio>
#include <cmath>

/* round a double */
double utils::nint(double d){
  return floor(d + 0.5);
}

/* modulus returning a positive number */
double utils::mod(double x, double y){
  double result = fmod(x,y);
  if(result<0) result += y;
  return result;
}

/* cross product of 2 vectors:
   c = a x b */
void utils::cross(double *a, double *b, double *c){
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/* transpose of a square matrix:
   AT = transpose(A) */
void utils::transpose(double *A, double *AT, int n){
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      AT[n*i+j] = A[n*j+i];
}

/* multiplication of 2 square matrices
   C = AB */
void utils::sqmatmul(double *A, double *B, double *C, int n){
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++){
      C[n*i+j] = 0.0;
      for(int k=0; k<n; k++)
	C[n*i+j] += A[n*i+k]*B[n*k+j];
    }
}

/* vector quaternion product */
void utils::vecquatprd(double *omega, double *q, double *qdot){
  qdot[0] = - omega[0]*q[1] - omega[1]*q[2] - omega[2]*q[3];
  qdot[1] =   omega[0]*q[0] + omega[1]*q[3] - omega[2]*q[2];
  qdot[2] =   omega[1]*q[0] - omega[0]*q[3] + omega[2]*q[1];
  qdot[3] =   omega[2]*q[0] + omega[0]*q[2] - omega[1]*q[1];
}

/* matrix-vector product with
   c = A*b, where b is 3x1, A is nx3 and c is nx1 */
void utils::matvecprd(double *A, double *b, double *c, int n){
  for(int i=0; i<n; i++){
    c[i] = 0.0;
    for(int j=0; j<3; j++)
      c[i] += A[3*i+j]*b[j];
  }
}

/* pretty-print square matrix */
void utils::printSqMatrix(double *A, int n){
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++)
      printf("%+e ", A[n*i+j]);
    printf("\n");
  }
}

/* pretty-print vector */
void utils::printVector(double *a, int n){
  printf("(%f", a[0]);
  for(int i=1; i<n; i++)
    printf(", %f", a[i]);
  printf(")\n");
}

/* create the '*' matrix out of a vector */
void utils::star(double *a, double *astar){
  astar[0] =  0.;
  astar[1] = -a[2];
  astar[2] =  a[1];
  astar[3] =  a[2];
  astar[4] =  0.;
  astar[5] = -a[0];
  astar[6] = -a[1];
  astar[7] =  a[0];
  astar[8] =  0.;
}
  
/* split string by whitespace */
using std::vector;
using std::string;
vector<string> utils::split(string &instr){
  string buf;
  vector<string> splitstr;
  std::stringstream ss(instr);
  while(ss>>buf)
    splitstr.push_back(buf);

  return splitstr;
}
