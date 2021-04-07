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
  int sizeA=800;
  int sizek=400;
  double deltaA=1./80.;
  double deltak=1./10.;

  
 public:
  MV(){  
    Phip_array = new double*[800];
    for (int i = 0; i < 800; i++) {
      Phip_array[i] = new double[400];
    }
  };

  ~MV(){  
    for(int i = 0; i < 800; ++i) {
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
