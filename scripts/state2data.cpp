#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>

enum{PERIODIC, EQUILIBRIUM, BOUNCEBACK};

int main(int argc, char* argv[]){
  // parse command line
  char *infile, *outfile;
  if(argc!= 5){
    std::cout<<"Usage: restart2data -f <infile> -o <outfile>\n";
    exit(0);
  }
  else{
    int ctr = 1;
    while(ctr<argc){
      if(ctr+1!=argc){
	if(strcmp(argv[ctr],"-f")==0){
	  infile = argv[ctr+1];
	  ctr += 2;
	}
	else if(strcmp(argv[ctr],"-o")==0){
	  outfile = argv[ctr+1];
	  ctr += 2;
	}
	else{
	  std::cerr<<"Not enough or invalid arguments\n";
	  std::cerr<<"Usage: restart2data -f <infile> -o <outfile>\n";
	  exit(0);
	}
      }
    }
  }

  // read input (binary) state file
  std::ifstream fin;
  fin.open(infile, std::ios::binary);
  if (!fin){
    std::cerr<<"ERROR: input file "<<infile<<" does not exist\n";
    exit(0);
  }

  fin.seekg(0, std::ios::end);
  unsigned int length = fin.tellg(); // length of file
  fin.seekg(0, std::ios::beg);

  double *buffer = new double[length/sizeof(double)];
  fin.read((char *)(buffer), length);
  fin.close();
  
  // write to output file
  std::ofstream fout;
  fout.open(outfile);
  int ctr = 0;

  // model velocity
  int nVelocity = static_cast<int>(buffer[ctr++]);

  // time step
  int tStep = static_cast<int>(buffer[ctr++]);
  fout<<tStep<<"\n";

  // domain dimensions
  int Lx = static_cast<int>(buffer[ctr++]);
  int Ly = static_cast<int>(buffer[ctr++]);
  int Lz = static_cast<int>(buffer[ctr++]);
  fout<<Lx<<" "<<Ly<<" "<<Lz<<"\n";

  // boundary type
  std::string bdryType[6];
  for(int i=0; i<6; i++){
    int tmp = static_cast<int>(buffer[ctr++]);
    if(tmp==PERIODIC)
      bdryType[i] = "p";
    else if(tmp==EQUILIBRIUM)
      bdryType[i] = "e";
    else if(tmp==BOUNCEBACK)
      bdryType[i] = "b";
    else{
      std::cerr<<"ERROR: domain bdry type "<<tmp<<" not recognized\n";
      exit(1);
    }
    fout<<bdryType[i]<<" ";
  }
  fout<<"\n";

  // boundary velocities
  double uE[3], uW[3], uN[3], uS[3], uU[3], uD[3];
  for(int i=0; i<3; i++){
    uE[i] = buffer[ctr++];
    fout<<uE[i]<<" ";
  }
  fout<<"\n";
  for(int i=0; i<3; i++){
    uW[i] = buffer[ctr++];
    fout<<uW[i]<<" ";
  }
  fout<<"\n";
  for(int i=0; i<3; i++){
    uN[i] = buffer[ctr++];
    fout<<uN[i]<<" ";
  }
  fout<<"\n";
  for(int i=0; i<3; i++){
    uS[i] = buffer[ctr++];
    fout<<uS[i]<<" ";
  }
  fout<<"\n";
  for(int i=0; i<3; i++){
    uU[i] = buffer[ctr++];
    fout<<uU[i]<<" ";
  }
  fout<<"\n";
  for(int i=0; i<3; i++){
    uD[i] = buffer[ctr++];
    fout<<uD[i]<<" ";
  }
  fout<<"\n";

  // LB parameters
  double nu, nu_b, extF[3];
  nu = buffer[ctr++];
  fout<<nu<<" ";
  nu_b = buffer[ctr++];
  fout<<nu_b<<" ";
  for(int i=0; i<3; i++){
    extF[i] = buffer[ctr++];
    fout<<extF[i]<<" ";
  }
  fout<<"\n";

  // MD parameters
  int nP; double mdTstep;
  nP = buffer[ctr++];
  mdTstep = buffer[ctr++];
  fout<<nP<<" ";
  fout<<mdTstep<<"\n";

  // finally, @ node: x, y, z, type, n(0,1,...,nVelocity-1)
  int x, y, z, type;
  double *n = new double[nVelocity];
  for(int xl=1; xl<=Lx; xl++){
    for(int yl=1; yl<=Ly; yl++){
      for(int zl=1; zl<=Lz; zl++){
	x = static_cast<int>(buffer[ctr++]);
	y = static_cast<int>(buffer[ctr++]);
	z = static_cast<int>(buffer[ctr++]);
	type = static_cast<int>(buffer[ctr++]);
	for(int k=0; k<nVelocity; k++)
	  n[k] = buffer[ctr++];

	fout<<x<<" "<<y<<" "<<z<<" "<<type<<" "<<;
	for(int k=0; k<nVelocity; k++)
	  fout<<n[k]<<" ";
	fout<<"\n";
      }
    }
  }
  fout<<"\n";

  fout.close();
  return 0;
}

