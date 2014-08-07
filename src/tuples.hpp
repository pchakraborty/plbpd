#include <iostream>

template<class T>
struct tuples3{
  T a, b, c;
  tuples3(){
  }
  tuples3(T u,T v,T w ){
    a = u; b = v; c = w;
  }
  tuples3(const tuples3& orig){
    a = orig.a; b = orig.b; c = orig.c;
  }
  tuples3& operator=(const tuples3& orig){
    a = orig.a; b = orig.b; c = orig.c;
    return *this; 
  }
};


template<class T>
std::ostream& operator<<(std::ostream& os, const tuples3<T>& t){
  // Save a triplet to os
  os << t.a << " " << t.b << " " << t.c;
  return os;
}
template<class T>
std::istream& operator>>(std::istream& is, tuples3<T>& t){
  // Load a tuples3 from is
  is >> t.a >> t.b >> t.c;
  return is;
}

template<class T>
struct tuples6{
  T a, b, c, d, e, f;
  tuples6(){
  }
  tuples6(T u,T v,T w,T x, T y, T z){
    a = u; b = v; c = w; d = x; e = y; f = z;
  }
  tuples6(const tuples6& orig){
    a = orig.a; b = orig.b; c = orig.c;
    d = orig.d; e = orig.e; f = orig.f;
  }
  tuples6& operator=(const tuples6& orig){
    a = orig.a; b = orig.b; c = orig.c;
    d = orig.d; e = orig.e; f = orig.f;
    return *this; 
  }
};


template<class T>
std::ostream& operator<<(std::ostream& os, const tuples6<T>& t){
  // Save a triplet to os
  os<<t.a<<" "<<t.b<<" "<<t.c<<" "<<t.d<<" "<<t.e<<" "<<t.f;
  return os;
}
template<class T>
std::istream& operator>>(std::istream& is, tuples6<T>& t){
  // Load a tuples3 from is
  is>>t.a>>t.b>>t.c>>t.d>>t.e>>t.f;
  return is;
}
