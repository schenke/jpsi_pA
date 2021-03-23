#ifndef TAINT_H_
#define TAINT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


using namespace std;


class TAInt{
 private:
  double TAgrid[200]; 
  double xgrid[200]; 
  void computeTAIntegral();
  static double rhoA(double z, void * params);

 public:
  TAInt(){
    computeTAIntegral();
      };
  ~TAInt(){};
  double returnTA(double R);

};

#endif // TAINT_H_
