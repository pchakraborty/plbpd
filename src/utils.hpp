#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <vector>
#include <sstream>

class utils{
public:
  // functions
  static double nint(double d);
  static double mod(double x, double y);
  static void cross(double *a, double *b, double *c);
  static void star(double *a, double *astar);
  static void transpose(double *A, double *AT, int n);
  static void sqmatmul(double *A, double *B, double *C, int n);
  static void matvecprd(double *A, double *b, double *c, int n);
  static void vecquatprd(double *omega, double *q, double *qdot);

  static void printSqMatrix(double *A, int n);
  static void printVector(double *a, int n);

  static std::vector<std::string> split(std::string &instr);
  
  template<class T>
  static inline std::string T2str(const T &t){
    std::stringstream ss;
    ss<<t;
    return ss.str();
  }

  template<class T>
  static inline T str2T(const std::string &str){
    std::istringstream buffer(str);
    T t;
    buffer>>t;
    return t;    
  }
};

#endif
