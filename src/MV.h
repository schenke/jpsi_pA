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
  static double BKintegrandForList(double z, void * params);
  double ** Phip_array;
  double *** Phip_arrayBK;
  int sizeA=3200;
  int sizek=1000;
  int sizey=18;
  double deltaA=1./80.;
  double deltak=1./10.;
  double deltay=1.;
  
 public:
  MV(){  
    Phip_array = new double*[sizeA];
    for (int i = 0; i < sizeA; i++) {
      Phip_array[i] = new double[sizek];
    }
    Phip_arrayBK = new double**[sizeA];
    for (int i = 0; i < sizeA; i++) {
      Phip_arrayBK[i] = new double*[sizek];
      for (int j = 0; j < sizek; j++) {
        Phip_arrayBK[i][j] = new double[sizey];
      }
    }
  };

  ~MV(){
    for(int i = 0; i < sizeA; ++i) {
      delete[] Phip_array[i];
    }
    delete[] Phip_array;
    for (int i = 0; i < sizeA; i++) {
      for (int j = 0; j < sizek; j++) {
        delete [] Phip_arrayBK[i][j];
      }
      delete [] Phip_arrayBK[i];
    }
    delete[] Phip_arrayBK;
  }; 
  
  void computePhip();
  void computePhipBK();
  double Phip(double k, double R, double Qs, double sizeFactor);
  double PhipFluc(double k, double Tp, double Qs, double sizeFactor);
  double PhipBKFluc(double k, double Tp, double Qs, double x);
  double Phit(double k, double TA, double Qs);
  double PhitBK(double k, double TA, double Qs,double x);
  double StF(double k, double TA, double Qs);
  double StFBK(double k, double TA, double Qs, double x);
  int writeTable();
  int writeTableBK();
  int writeTableText();
  int readTable();
  int readTableBK();

};

#endif // MV_H_
