#ifndef ERROR_HPP
#define ERROR_HPP

#include "pointers.hpp"

class Error: protected Pointers{
public:
  Error(class PLBPD *);
  ~Error();
  
  void all(string &);
  void warning(string &);
  void notify(string &);
};

#endif
