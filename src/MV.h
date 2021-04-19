#ifndef MV_H_
#define MV_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;


class MV{
 private:
  static double MVintegrandForList(double z, void * params);
  double ** Phip_array;
  int sizeA=1600;
  int sizek=1000;
  double deltaA=1./80.;
  double deltak=1./10.;
  
 public:
  MV(){  
    Phip_array = new double*[sizeA];
    for (int i = 0; i < sizeA; i++) {
      Phip_array[i] = new double[sizek];
    }
  };

  ~MV(){
    for(int i = 0; i < sizeA; ++i) {
      delete[] Phip_array[i];
    }
    delete[] Phip_array;
  }; 
  
  void computePhip();
  double Phip(double k, double R, double Qs);
  double Phit(double k, double TA, double Qs);
  double StF(double k, double TA, double Qs);
  int writeTable();
  int readTable();

};

#endif // MV_H_
