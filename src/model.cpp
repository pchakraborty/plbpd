#include "model.hpp"
#include "D3Q19.hpp"
#include "D3Q27.hpp"
#include "input.hpp"
#include "error.hpp"

/******************************************************************************/

Model::Model(PLBPD *plbpd): Pointers(plbpd){
  fluid = NULL;
  thermal = NULL;
}

/******************************************************************************/

Model::~Model(){
  if(fluid) delete fluid;
  if(thermal) delete thermal;
}

/******************************************************************************/

void Model::setup(){
//   // 2D
//   if(input->Ly==1){
//     fluid = new Var;
//     fluid->nVelocity = d2q9_nVelocity;
//     fluid->cs2 = d2q9_cs2;
//     fluid->c = d2q9_c;
//     fluid->w = d2q9_w;
//     fluid->n_unk = d2q9_n_unk;
//     fluid->unkE = d2q9_unk_E;
//     fluid->unkW = d2q9_unk_W;
//     fluid->unkN = d2q9_unk_N;
//     fluid->unkS = d2q9_unk_S;
//     fluid->unkU = d2q9_unk_U;
//     fluid->unkD = d2q9_unk_D;
//     fluid->reverse = d2q9_reverse;
//     //mirror = d2q9_mirror;
    
//     if(input->heat_xfer){
//       thermal = new Var;
//       thermal->nVelocity = d2q9_nVelocity;
//       thermal->cs2 = d2q9_cs2;
//       thermal->c = d2q9_c;
//       thermal->w = d2q9_w;
//       //coeff = d2q9_coeffs;
//       thermal->n_unk = d2q9_n_unk;
//       thermal->unkE = d2q9_unk_E;
//       thermal->unkW = d2q9_unk_W;
//       thermal->unkN = d2q9_unk_N;
//       thermal->unkS = d2q9_unk_S;
//       thermal->unkU = d2q9_unk_U;
//       thermal->unkD = d2q9_unk_D;
//       thermal->reverse = d2q9_reverse;
//       //mirror = d2q9_mirror;
//     }
//   }

//   // 3D
//   else if(input->Ly>1){
  if(input->Ly>=1){ // 2D or 3D
    fluid = new Var;
    fluid->nVelocity = d3q19_nVelocity;
    fluid->cs2 = d3q19_cs2;
    fluid->c = d3q19_c;
    fluid->w = d3q19_w;
    //coeff = d3q19_coeffs;
    fluid->n_unk = d3q19_n_unk;
    fluid->unkE = d3q19_unk_E;
    fluid->unkW = d3q19_unk_W;
    fluid->unkN = d3q19_unk_N;
    fluid->unkS = d3q19_unk_S;
    fluid->unkU = d3q19_unk_U;
    fluid->unkD = d3q19_unk_D;
    fluid->reverse = d3q19_reverse;
    //mirror = d3q19_mirror;
    
    if(input->heat_xfer){
      thermal = new Var;
      thermal->nVelocity = d3q27_nVelocity;
      thermal->cs2 = d3q27_cs2;
      thermal->c = d3q27_c;
      thermal->w = d3q27_w;
      //coeff = d3q27_coeffs;
      thermal->n_unk = d3q27_n_unk;
      thermal->unkE = d3q27_unk_E;
      thermal->unkW = d3q27_unk_W;
      thermal->unkN = d3q27_unk_N;
      thermal->unkS = d3q27_unk_S;
      thermal->unkU = d3q27_unk_U;
      thermal->unkD = d3q27_unk_D;
      thermal->reverse = d3q27_reverse;
      //mirror = d3q27_mirror;
    }
  }
  else{ // input->Ly<1
    string str = "domain dimension in the y-direction can either be";
    str += " =1 (2D) or >1 (3D)";
    error->all(str);
  }
}

/******************************************************************************/
