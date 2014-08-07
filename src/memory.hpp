#ifndef MEMORY_HPP
#define MEMORY_HPP

#include "pointers.hpp"
#include "error.hpp"

class Memory: protected Pointers{
public:
  Memory(PLBPD *plbpd): Pointers(plbpd){}
  ~Memory(){}

  template<class T>
  T ****create4Darray(int n1, int n2, int n3, int n4){

    T *data = NULL, **cube = NULL, ***plane = NULL, ****array = NULL;

    try{
      data = new T[n1*n2*n3*n4];
      cube = new T*[n1*n2*n3];
      plane = new T**[n1*n2];
      array = new T***[n1];
    }
    catch(std::bad_alloc &xa){
      string str("Memory::create4Darray() failed to allocate memory");
      error->all(str);
    }

    int n = 0;
    for(int i=0; i<n1; i++){
      array[i] = &plane[i*n2];
      for(int j=0; j<n2; j++){
	plane[i*n2+j] = &cube[i*n2*n3+j*n3];
	for(int k=0; k<n3; k++){
	  cube[i*n2*n3+j*n3+k] = &data[n];
	  n += n4;
	}
      }
    }
    
    return array;
  }

  template<class T>  
  void delete4Darray(T ****array){
    if(array==NULL) return;
    delete [] array[0][0][0];
    delete [] array[0][0];
    delete [] array[0];
    delete [] array;
  }

  template<class T>
  T ***create3Darray(int n1, int n2, int n3){

    T *data = NULL, **plane = NULL, ***array = NULL;
    try{
      data = new T[n1*n2*n3];
      plane = new T*[n1*n2];
      array = new T**[n1];
    }
    catch(std::bad_alloc &xa){
      string str("Memory::create3Darray() failed to allocate memory");
      error->all(str);
    }
    
    int n = 0;
    for(int i=0; i<n1; i++){
      array[i] = &plane[i*n2];
      for(int j=0; j<n2; j++){
	plane[i*n2+j] = &data[n];
	n += n3;
      }
    }
    
    return array;
  }

  template<class T>
  void delete3Darray(T ***array){
    if(array==NULL) return;
    delete [] array[0][0];
    delete [] array[0];
    delete [] array;
  }
  
  template<class T>
  T **create2Darray(int nrows, int ncols){

    T *data = NULL, **array = NULL;
    try{
      data = new T[nrows*ncols];
      array = new T*[nrows];
    }
    catch(std::bad_alloc &xa){
      string str("Memory::create2Darray() failed to allocate memory");
      error->all(str);
    }

    int n = 0;
    for(int i=0; i<nrows; i++){
      array[i] = &data[n];
      n += ncols;
    }
    
    return array;
  }

  template<class T>
  void delete2Darray(T **array){
    if(array==NULL) return;
    delete [] array[0];
    delete [] array;
  }    
};

#endif
