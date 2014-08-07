// data structure for primaries

/*
  CURRENTLY, OTHER THAN THE ID, THIS CLASS HAS 7 DOUBLES,
  r, c[3] and cb[3]. IF EVER THIS CHANGES, YOU NEED TO EDIT
  THE VARIABLES 'NPRMRYDBLES' IN 'restart.hpp'. ALSO, THE
  'write agglomerate data' SECTION HAS TO BE MODIFIED ACCORDINGLY
*/

#ifndef PRIMARY_HPP
#define PRIMARY_HPP

class Primary{
public:
  unsigned int ID;       // unique ID
  double r;              // primary (sphere) radius
  double c[3];           // primary centre (global space)
  double cb[3];          // primary centre (body space) [DOES NOT CHANGE]

  Primary(){
    r = -1.0;
    for(int i=0; i<3; i++){
      c[i] = -1.0;
      cb[i] = -1.0;
    }
  }
};

#endif
