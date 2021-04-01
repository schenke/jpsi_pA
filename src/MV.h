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
  
 public:
  MV(){  
    Phip_array = new double*[200];
    for (int i = 0; i < 200; i++) {
      Phip_array[i] = new double[200];
    }
  };

  ~MV(){  
    for(int i = 0; i < 200; ++i) {
      delete[] Phip_array[i];  
    }
    delete[] Phip_array;  
  };
  
  void computePhip();
  double Phip(double k, double R, double Qs);

};

#endif // MV_H_
