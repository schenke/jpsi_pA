// Parameters.h is part of the simpleMCGlauber model.
// Copyright (C) 2014 Bjoern Schenke.

#ifndef Parameters_H
#define Parameters_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "mpi.h"


using namespace std;

class Parameters 
{
 private:
  //constants
  double myhbarc;

  // switches:
  int gaussianWounding;
  int timeForSeed;
  int useQuarks;
  int useEnergyDependentCrossSection;

  //integers:
  int MPIrank;
  int MPIsize;
  long seed;
  int outputNumberOfTransverseCells;

  //doubles:
  double b;
  double sigmaNN;
  double sigmaQQ;
  double sigmaNNinel;
  double sigmaQQinel;
  double roots; // center of mass energy
  double outputMaximalTransverseSize;
  double transverseGaussianSmearingWidth;

  //strings:
  string projectile;
  string target;

  string temp;

 public:

  // constructor:
  Parameters() 
    {
      sethbarc(0.1973269631);   //hbar c in GeV*fm
    }

  // functions to access the private variables:
  string strWord(int index, string line);

  //strings
  void setProjectile(string x){projectile=x;}
  string getProjectile(){return projectile;}
  void setTarget(string x){target=x;}
  string getTarget(){return target;}

  //bools
  bool readParameter(string parameter);
  void setParameters(string Target);

  //longs
  void setSeed(long x){seed=x;}
  long getSeed(){return seed;}
  
  //ints
  void setMPIRank(int x){MPIrank=x;}
  int getMPIRank(){return MPIrank;}
  void setMPISize(int x){MPIsize=x;}
  int getMPISize(){return MPIsize;}
  void setGaussianWounding(int x){gaussianWounding=x;}
  int getGaussianWounding(){return gaussianWounding;}
  void setTimeForSeed(int x){timeForSeed=x;}
  int getTimeForSeed(){return timeForSeed;}
  void setUseQuarks(int x){useQuarks=x;}
  int getUseQuarks(){return useQuarks;}
  void setUseEnergyDependentCrossSection(int x){useEnergyDependentCrossSection=x;}
  int getUseEnergyDependentCrossSection(){return useEnergyDependentCrossSection;}
  void setOutputNumberOfTransverseCells(int x){outputNumberOfTransverseCells=x;}
  int getOutputNumberOfTransverseCells(){return outputNumberOfTransverseCells;}


  //doubles
  void sethbarc(double x){myhbarc=x;}
  double gethbarc(){return myhbarc;}
  void setb(double x){b=x;}
  double getb(){return b;}
  void setSigmaNN(double x){sigmaNN=x;}
  double getSigmaNN(){return sigmaNN;}
  void setSigmaQQ(double x){sigmaQQ=x;}
  double getSigmaQQ(){return sigmaQQ;}
  void setSigmaNNinel(double x){sigmaNNinel=x;}
  double getSigmaNNinel(){return sigmaNNinel;}
  void setSigmaQQinel(double x){sigmaQQinel=x;}
  double getSigmaQQinel(){return sigmaQQinel;}
  void setRoots(double x){roots=x;}
  double getRoots(){return roots;}
  void setOutputMaximalTransverseSize(double x){outputMaximalTransverseSize=x;}
  double getOutputMaximalTransverseSize(){return outputMaximalTransverseSize;}
  void setTransverseGaussianSmearingWidth(double x){transverseGaussianSmearingWidth=x;}
  double getTransverseGaussianSmearingWidth(){return transverseGaussianSmearingWidth;}



};
#endif // Parameters_H
