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
  Glauber(Parameters *inParam);
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
  int collide(Random *random);
  int getNpart();
  double FermiDistribution(Nucleus *nucleus, double r);
  void init(Random *random);
  void makeHistograms(Random *random);
  Nucleon sampleRho(Nucleus *nucleus, Random *random); 
  Quark sampleQuark(Random *random);
  void outputNucleonPositions();
  void outputQuarkPositions();
  double ExponentialDistribution(double a, double r);
  void makeNuclei(Random *random);
  int getAverageNumberOfCollisions();
  vector<complex <double> > eccentricity(int n, string folder);
  void standard_eccentricity(int event, string folder);

  double unitStep(double a);
  
  Nucleus getTarget(){return Target;}
  Nucleus getProjectile(){return Projectile;}

  void computeSigmaNN(Random* random);


 private:
  Parameters *param;
  Nucleus Target;
  Nucleus Projectile;

  double *hadronBins;

  gsl_rng * gslRan;

  int numberOfQuarks;

 };
#endif

