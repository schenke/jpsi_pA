#ifndef Glauber_h //avoids multiple inclusions of the header file
#define Glauber_h

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "Random.h"
#include "Parameters.h"
#include <algorithm> // for sort
#include "mpi.h"
#include "gsl/gsl_rng.h"

using namespace std;

struct Quark
{
  double x;
  double y;
  double z;
  int collided;
  double rapidity;
  vector<int> collidedWith;
  vector<int> collidedWithQuark;
};

struct Nucleon 
{
  double x;
  double y;
  double z;
  int collided;
  double rapidity;
  vector<Quark> quarkList;
  vector<int> collidedWith;
};

class Glauber{

 public:
  Glauber(Parameters *inParam, double inwidth);
  ~Glauber();

  typedef double (*ptr_func)(double);
  
  typedef struct nucleus 
  {
    string name;
    double A;
    double Z;
    int DensityFunc;
    double w_WS;
    double a_WS;
    double R_WS;
    double rho_WS;
    vector<Nucleon> nucleonList;
   
  } Nucleus;

  int setNucleusParameters(Nucleus *nucleus, string name);
  double FermiDistribution(Nucleus *nucleus, double r);
  void init(Random *random);
  void makeHistograms(Random *random);
  Nucleon sampleRho(Nucleus *nucleus, Random *random); 
  Quark sampleQuark(Random *random);
  double ExponentialDistribution(double a, double r);
  void makeNuclei(Random *random, double Bp, double Bq);
  void generateProtonTp(Nucleus *nuc, Random *random, double Bp, double Bq);
  void generateNucleusTA(Nucleus *nuc, Random *random, double Bp, double Bq);
  double returnNucleusTA(double x, double y);
  double returnProtonTp(double x, double y);

  int getA();
  Nucleus getTarget(){return Target;}
  Nucleus getProjectile(){return Projectile;}
  //  double returnTA2D(double x, double y);


 private:
  Parameters *param;
  const double PI = 3.14159265358979323846;
  Nucleus Target;
  Nucleus Projectile;

  double *hadronBins;

  gsl_rng * gslRan;

  int numberOfQuarks;
  double width;
  double TAgrid2D[200][200]; 
  double Tpgrid2D[40][40]; 
  
 };
#endif

