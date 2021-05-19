// Jpsi and gluon cross section calculation in MV or GBW models 
// including fluctuating geometries 
// 2021 Bjoern Schenke

#include "mpi.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
//#include <vector>
//#include <ctype.h>

#include <gsl/gsl_integration.h>
//#include <gsl/gsl_interp.h>
//#include <gsl/gsl_spline.h>
#include "cuba.h"
#include "pretty_ostream.h"
#include "hard.h"
#include "TAInt.h"
#include "MV.h"
#include "Glauber.h"
#include "nrqcd.h"

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

void display_logo();

// Constants to be used
// MV.cpp has constants too, make sure they are the same... needs to be improved
namespace constants {
  const double PI = 3.14159265358979323846;
  const int Nc = 3;
  const double hbarc = 0.1973269804;
  const double CA = double(Nc);
  const double CF = (double(Nc)*double(Nc) - 1.)/(2.*double(Nc));
  const double alphas = 0.3;
  const double Bp = 4.; // == R_p = 0.4 fm
  const double Bt = 1061; // == R_t = 1.1 * A^{1/3} fm ~6.5 fm
  const double mD = 1.864;
  const double mc = 1.275; //vary? 1.4?
  const double mJPsi = 3.096916;
  const double x0 = 0.000041;
  //const double lambdaSpeed = 0.277;
  //const double x0 = 0.00005;
  const double lambdaSpeedp = 0.277;
  const double lambdaSpeedA = 0.277;
  const double prefactor = 0.56;
  //const double prefactor = 0.7;
  const double roots = 8160.;
  const double ldme_singlet = 1.16/2./double(Nc); // GeV^3
  const double ldme_octet_s10 = 0.089; // +- 0.0098 GeV^3
  const double ldme_octet_s13 = 0.0030; // +- 0.00012 GeV^3
  const double ldme_octet_p3j = 0.0056*mc*mc; // (+- 0.0021 GeV^3) [GeV^5]

  const double bdep_p = 1.87; // Eq. 114 notes. We plug it in the MV.cpp
  const double bdep_A = 4.13; // Eq. 118 notes  
  const double bindep_A = 2.77; // Eqs. 118 and 123 notes == 4.13*0.67
  const double bdep_fluc_p = 1.87; // Eq. 114 notes
  const double bdep_fluc_A = 1.87; // same as for the proton
}

// Parameters that need to be passed to the integrand
struct params {
  double pe;
  double k;
  double Qs;
  double Qsp;
  double QsA;
  double lambda;
  double Y;
  double m;
  double PT;
  double bx;
  double by;
  double yq;
  double yp;
  TAInt *TAclass;
  MV *mv;
  Glauber *glauberClass; 
  double protonSizeFactor;
  int BK;
};

struct kinqq {
  double p;
  double phip;
  double q;
  double phiq;
  double yp;
  double yq;
};

struct kinPair {
  double M;
  double PT;
  double Y;
  double qtilde;
  double phi;
  double m;
};

// Conversion of kinematic variables
kinqq convert(kinPair input) {

  kinqq output;

  double M = input.M;
  double PT = input.PT;
  double Y = input.Y;
  double qtilde = input.qtilde;
  double phi = input.phi;
  double m = input.m;
  
  double betax = PT/sqrt(M*M+PT*PT);
  double gammax = 1./sqrt(1.-betax*betax);
  double betaz = tanh(Y);
  double gammaz = 1./sqrt(1.-betaz*betaz);

  double pE = gammaz*(gammax*(M/2.-betax*qtilde*cos(phi))
                      -betaz*sqrt(M*M/4.-m*m-qtilde*qtilde));
  double px = gammax*(betax*M/2.-qtilde*cos(phi));
  double py = -qtilde*sin(phi);
  double pz = gammaz*(betaz*gammax*(M/2.-betax*qtilde*cos(phi))
                      -sqrt(M*M/4.-m*m-qtilde*qtilde));
  
  double qE = gammaz*(gammax*(M/2.+betax*qtilde*cos(phi))
                      +betaz*sqrt(M*M/4.-m*m-qtilde*qtilde));
  double qx = gammax*(betax*M/2.+qtilde*cos(phi));
  double qy = qtilde*sin(phi);
  double qz = gammaz*(betaz*gammax*(M/2.+betax*qtilde*cos(phi))
                      +sqrt(M*M/4.-m*m-qtilde*qtilde));

  output.p = sqrt(px*px+py*py);
  output.phip = atan2(py,px);
  output.q = sqrt(qx*qx+qy*qy);
  output.phiq = atan2(qy,qx);
  output.yp = 0.5*log((pE+pz)/(pE-pz));
  output.yq = 0.5*log((qE+qz)/(qE-qz));
  
  return output;
}
  
double returnTA2D(double x, double y, Glauber *glauberClass){
  return glauberClass->returnNucleusTA(x, y);
}

double returnTp2D(double x, double y, Glauber *glauberClass){
  return glauberClass->returnProtonTp(x, y);
}

double returnTA(double R, TAInt *TAclass){
return TAclass->returnTA(R);
}

// Unintegrated gluon distribution for the proton in GBW
double PhipGBW(double k, double R, double Qs){
  return constants::CF*k*k*constants::Nc*constants::PI/constants::alphas/constants::CA/Qs/Qs*exp(R*R/(2.*constants::Bp))
    *exp(-constants::CF*exp(R*R/(2.*constants::Bp))*k*k/(constants::CA*Qs*Qs));
}

// choose between MV and GBW - should make this choice a parameter of course
double PhipFluc(double k, double Tp, double Qs, double sizeFactor, MV *mv, int BK, double x){
  if(BK){
    return mv->PhipBKFluc(k, Tp, x);
  }
  else{
    return mv->PhipFluc(k, Tp, Qs, sizeFactor);
  }
  //  return mv->Phip(k, RorTp, Qs, sizeFactor);
  //return PhipGBW(k, R, Qs);
}

double Phip(double k, double R, double Qs, double sizeFactor, MV *mv, int BK, double x){
  if(BK){
    return mv->PhipBK(k, R, sizeFactor,x);
  }
  else{
    return mv->Phip(k, R, Qs, sizeFactor);
  }
}

// Unintegrated gluon distribution for lead
double PhitGBW(double k, double TA, double Qs){
  return constants::PI*k*k*constants::Nc/constants::alphas*constants::CF*exp(-constants::CF*k*k/(constants::CA*Qs*Qs*TA))/(constants::CA*Qs*Qs*TA);
}

// choose between MV and GBW - should make this choice a parameter of course
double Phit(double k, double TA, double Qs, MV *mv, int BK, double x){
  if(BK){
    return mv->PhitBK(k, TA, x);
  }
  else{
    return mv->Phit(k, TA, Qs);
  }
  //return PhitGBW(k, TA, Qs);
}

// FT of fundamental S of the target
double StFGBW(double k, double TA, double Qs){
  return 4.*constants::PI*exp(-(k*k/(Qs*Qs*TA)))/(Qs*Qs*TA);
} 

// choose between MV and GBW - should make this choice a parameter of course
double StF(double k, double TA, double Qs, MV *mv, int BK, double x){
  if(BK){
    return mv->StFBK(k, TA, x);
  }
  else{
    return mv->StF(k, TA, Qs);
  }
//return StFGBW(k, TA, Qs);
}

static int JPsiIntegrandNRQCDCsNob(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobcs4k xx[0]
#define nobcs4phik xx[1]
#define nobcs4k1 xx[2]
#define nobcs4phik1 xx[3]
#define nobcs4kprime xx[4]
#define nobcs4phikprime xx[5]
#define nobcs4p xx[6]
#define nobcs4phip xx[7]
#define f ff[0]

  double kscale = 30.;
  double pscale = 30.;

  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double p = nobcs4p*pscale;
  double phip = nobcs4phip*2.*constants::PI;
  double k = nobcs4k*kscale;
  double phik = nobcs4phik*2.*constants::PI;
  double k1 = nobcs4k1*kscale;
  double phik1 = nobcs4phik1*2.*constants::PI;
  double kprime = nobcs4kprime*kscale;
  double phikprime = nobcs4phikprime*2.*constants::PI;
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/constants::roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double kprimex = kprime*cos(phikprime);
  double kprimey = kprime*sin(phikprime);
  double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
  double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
  double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                    +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
  double phi_pminuskminusk1minuskprime = atan2(pminuskminusk1minuskprimey,pminuskminusk1minuskprimex);

  
  double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);
 
  double myTA = 1.*constants::bindep_A;
 
  if(pminuskminusk1minuskprime>30.){
    f=0.;
  }
  else{
    f = constants::alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA)*factorxA*StF(kprime,myTA,QsA,mv,BK,xA)*factorxA*StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA)*factorxA)
      *2.*constants::PI*constants::Bp
      *2.*constants::PI*constants::Bt
      *p*pscale*2.*constants::PI
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
    // d2R
    // d2b
    // d2p
    // d2k
    // d2k1
    // d2kprime
  }
  }
  return 0;
} 

// Integrand for the combined J/Psi integral in the color singlet channel
static int JPsiIntegrandNRQCDCs(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define csR xx[0]
#define csphiR xx[1]
#define csb xx[2]
#define csphib xx[3]
#define cs4k xx[4]
#define cs4phik xx[5]
#define cs4k1 xx[6]
#define cs4phik1 xx[7]
#define cs4kprime xx[8]
#define cs4phikprime xx[9]
#define cs4p xx[10]
#define cs4phip xx[11]


  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!

  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double R = csR*Rscale;
  double b = csb*bscale;
  double p = cs4p*pscale;
  double phip = cs4phip*2.*constants::PI;
  double k = cs4k*kscale;
  double phik = cs4phik*2.*constants::PI;
  double k1 = cs4k1*kscale;
  double phik1 = cs4phik1*2.*constants::PI;
  double kprime = cs4kprime*kscale;
  double phikprime = cs4phikprime*2.*constants::PI;
  double phiR = csphiR*2*constants::PI;
  double phib = csphib*2*constants::PI;
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/constants::roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double kprimex = kprime*cos(phikprime);
  double kprimey = kprime*sin(phikprime);
  double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
  double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
  double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                    +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
  double phi_pminuskminusk1minuskprime = atan2(pminuskminusk1minuskprimey,pminuskminusk1minuskprimex);

  double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
  
  double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);
 
  double myTA = returnTA(Rminusb,TAclass)*constants::bdep_A; //(2.37=0.4*(208)^(1/3))
  
  if(pminuskminusk1minuskprime>30.){
    f=0.;
  }
  else{
    f = constants::alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA)*factorxA*StF(kprime,myTA,QsA,mv,BK,xA)*factorxA*StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA)*factorxA)
      *R*Rscale*2.*constants::PI
      *b*bscale*2.*constants::PI
      *p*pscale*2.*constants::PI
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
    // d2R
    // d2b
    // d2p
    // d2k
    // d2k1
    // d2kprime
   }
  }
  return 0;
} 

static int JPsiIntegrandNRQCDCsNoPT(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define csnoptR xx[0]
#define csnoptphiR xx[1]
#define csnoptb xx[2]
#define csnoptphib xx[3]
#define csnopt4k xx[4]
#define csnopt4phik xx[5]
#define csnopt4k1 xx[6]
#define csnopt4phik1 xx[7]
#define csnopt4kprime xx[8]
#define csnopt4phikprime xx[9]
#define csnopt4phip xx[10]

  double kscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!

  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double R = csnoptR*Rscale;
  double b = csnoptb*bscale;
  double p = PT;
  double phip = csnopt4phip*2.*constants::PI;
  double k = csnopt4k*kscale;
  double phik = csnopt4phik*2.*constants::PI;
  double k1 = csnopt4k1*kscale;
  double phik1 = csnopt4phik1*2.*constants::PI;
  double kprime = csnopt4kprime*kscale;
  double phikprime = csnopt4phikprime*2.*constants::PI;
  double phiR = csnoptphiR*2*constants::PI;
  double phib = csnoptphib*2*constants::PI;
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/constants::roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double kprimex = kprime*cos(phikprime);
  double kprimey = kprime*sin(phikprime);
  double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
  double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
  double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                    +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
  double phi_pminuskminusk1minuskprime = atan2(pminuskminusk1minuskprimey,pminuskminusk1minuskprimex);

  double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
  
  double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);
 
  double myTA = returnTA(Rminusb,TAclass)*constants::bdep_A; //(2.37=0.4*(208)^(1/3))
  
  if(pminuskminusk1minuskprime>30.){
    f=0.;
  }
  else{
    f = constants::alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA)*factorxA*StF(kprime,myTA,QsA,mv,BK,xA)*factorxA*StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA)*factorxA)
      *R*Rscale*2.*constants::PI
      *b*bscale*2.*constants::PI
      *p*2.*constants::PI
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
    // d2R
    // d2b
    // ptdphip
    // d2k
    // d2k1
    // d2kprime
   }
  }
  return 0;
} 

// Integrand for the combined J/Psi integral in the color singlet channel with fluctuations
static int JPsiIntegrandNRQCDCsFluc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define fcsRx xx[0]
#define fcsRy xx[1]
#define fcs4k xx[2]
#define fcs4phik xx[3]
#define fcs4k1 xx[4]
#define fcs4phik1 xx[5]
#define fcs4kprime xx[6]
#define fcs4phikprime xx[7]
#define fcs4p xx[8]
#define fcs4phip xx[9]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)

  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double Rx = fcsRx*Rscale-Rscale/2.;
  double Ry = fcsRy*Rscale-Rscale/2.;
  double p = fcs4p*pscale;
  double phip = fcs4phip*2.*constants::PI;
  double k = fcs4k*kscale;
  double phik = fcs4phik*2.*constants::PI;
  double k1 = fcs4k1*kscale;
  double phik1 = fcs4phik1*2.*constants::PI;
  double kprime = fcs4kprime*kscale;
  double phikprime = fcs4phikprime*2.*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/constants::roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double kprimex = kprime*cos(phikprime);
  double kprimey = kprime*sin(phikprime);

  double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
  double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
  double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                    +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
  
  double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);
  
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass)*constants::bdep_fluc_A;
  double Tp = returnTp2D(Rx,Ry,glauberClass)*constants::bdep_fluc_p;
  
  if(pminuskminusk1minuskprime>30.){
    f=0.;
  }
  else
    {
      f = constants::alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp)*factorxp/(k1*k1)*H_cs
        *(StF(k,TA,QsA,mv, BK, xA)*factorxA*StF(kprime,TA,QsA,mv, BK, xA)*factorxA*StF(pminuskminusk1minuskprime,TA,QsA,mv, BK, xA)*factorxA)
        *Rscale*Rscale
        *p*pscale*2.*constants::PI
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI
        *kprime*kscale*2.*constants::PI; 
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // dRxdRy
  // d2p
  // d2k
  // d2k1
  // d2kprime
    }
  }
  return 0;
} 

static int JPsiIntegrandNRQCDCoNob(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobco4k xx[0]
#define nobco4phik xx[1]
#define nobco4k1 xx[2]
#define nobco4phik1 xx[3]
#define nobco4p xx[4]
#define nobco4phip xx[5]

  double kscale = 30.;
  double pscale = 30.;

  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double p = nobco4p*pscale;
  double phip = nobco4phip*2.*constants::PI;
  double k = nobco4k*kscale;
  double phik = nobco4phik*2.*constants::PI;
  double k1 = nobco4k1*kscale;
  double phik1 = nobco4phik1*2.*constants::PI;
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/constants::roots;
  
  //!!
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pminuskminusk1x = px-kx-k1x;
  double pminuskminusk1y = py-ky-k1y;
  double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  double phi_pminuskminusk1 = atan2(pminuskminusk1y,pminuskminusk1x);
  
  double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
  double myTA = 1.*constants::bindep_A;
  
  f = constants::alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp)*factorxp/(k1*k1)*H_co
    *(StF(k,myTA,QsA,mv, BK, xA)*factorxA*StF(pminuskminusk1,myTA,QsA,mv, BK, xA)*factorxA)
    *2.*constants::PI*constants::Bp
    *2.*constants::PI*constants::Bt
    *p*pscale*2.*constants::PI
    *k*kscale*2.*constants::PI
    *k1*kscale*2.*constants::PI; 
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // d2p
  // d2k
  // d2k1
   
  }
  return 0;
}

// Integrand for the combined J/Psi integral in the color octet channel
static int JPsiIntegrandNRQCDCo(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define coR xx[0]
#define cophiR xx[1]
#define cob xx[2]
#define cophib xx[3]
#define co4k xx[4]
#define co4phik xx[5]
#define co4k1 xx[6]
#define co4phik1 xx[7]
#define co4p xx[8]
#define co4phip xx[9]

  double kscale = 15.;
  double pscale = 15.;
  double Rscale = 2./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!
  
  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double R = coR*Rscale;
  double b = cob*bscale;
  double p = co4p*pscale;
  double phip = co4phip*2.*constants::PI;
  double k = co4k*kscale;
  double phik = co4phik*2.*constants::PI;
  double k1 = co4k1*kscale;
  double phik1 = co4phik1*2.*constants::PI;
  double phiR = cophiR*2*constants::PI;
  double phib = cophib*2*constants::PI;
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/constants::roots;
  
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pminuskminusk1x = px-kx-k1x;
  double pminuskminusk1y = py-ky-k1y;
  double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  double phi_pminuskminusk1 = atan2(pminuskminusk1y,pminuskminusk1x);

  double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
  
  double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
  double myTA = returnTA(Rminusb,TAclass)*constants::bdep_A; //(2.37=0.4*(208)^(1/3));
  
  f = constants::alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp)*factorxp/(k1*k1)*H_co
    *(StF(k,myTA,QsA,mv, BK, xA)*factorxA*StF(pminuskminusk1,myTA,QsA,mv, BK, xA)*factorxA)
    *R*Rscale*2.*constants::PI
    *b*bscale*2.*constants::PI
    *p*pscale*2.*constants::PI
    *k*kscale*2.*constants::PI
    *k1*kscale*2.*constants::PI; 
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // d2p
  // d2k
  // d2k1
   
  }
  return 0;
} 

// Integrand for the combined J/Psi integral in the color octet channel
static int JPsiIntegrandNRQCDCoNoPT(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define conoptR xx[0]
#define conoptphiR xx[1]
#define conoptb xx[2]
#define conoptphib xx[3]
#define conopt4k xx[4]
#define conopt4phik xx[5]
#define conopt4k1 xx[6]
#define conopt4phik1 xx[7]
#define conopt4phip xx[8]

  double kscale = 30.;
  double Rscale = 2./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!
  
  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double R = conoptR*Rscale;
  double b = conoptb*bscale;
  double p = PT;
  double phip = conopt4phip*2.*constants::PI;
  double k = conopt4k*kscale;
  double phik = conopt4phik*2.*constants::PI;
  double k1 = conopt4k1*kscale;
  double phik1 = conopt4phik1*2.*constants::PI;
  double phiR = conoptphiR*2*constants::PI;
  double phib = conoptphib*2*constants::PI;
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/constants::roots;
  
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pminuskminusk1x = px-kx-k1x;
  double pminuskminusk1y = py-ky-k1y;
  double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  double phi_pminuskminusk1 = atan2(pminuskminusk1y,pminuskminusk1x);

  double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
  
  double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
  double myTA = returnTA(Rminusb,TAclass)*constants::bdep_A; //(2.37=0.4*(208)^(1/3));
  
  f = constants::alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp)*factorxp/(k1*k1)*H_co
    *(StF(k,myTA,QsA,mv, BK, xA)*factorxA*StF(pminuskminusk1,myTA,QsA,mv, BK, xA)*factorxA)
    *R*Rscale*2.*constants::PI
    *b*bscale*2.*constants::PI
    *p*2.*constants::PI
    *k*kscale*2.*constants::PI
    *k1*kscale*2.*constants::PI; 
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // dp_phi
  // d2k
  // d2k1
   
  }
  return 0;
} 

// Integrand for the combined J/Psi integral in the color octet channel for fluctuations
static int JPsiIntegrandNRQCDCoFluc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define fcoRx xx[0]
#define fcoRy xx[1]
#define fco4k xx[2]
#define fco4phik xx[3]
#define fco4k1 xx[4]
#define fco4phik1 xx[5]
#define fco4p xx[6]
#define fco4phip xx[7]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)

  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double Rx = fcoRx*Rscale-Rscale/2.;
  double Ry = fcoRy*Rscale-Rscale/2.;
  double p = fco4p*pscale;
  double phip = fco4phip*2.*constants::PI;
  double k = fco4k*kscale;
  double phik = fco4phik*2.*constants::PI;
  double k1 = fco4k1*kscale;
  double phik1 = fco4phik1*2.*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/constants::roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/constants::roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pminuskminusk1x = px-kx-k1x;
  double pminuskminusk1y = py-ky-k1y;
  double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);

  double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass)*constants::bdep_fluc_A;
  double Tp = returnTp2D(Rx,Ry,glauberClass)*constants::bdep_fluc_p;
 
  if(pminuskminusk1>30.){
    f=0.;
  }
  else
    {
      f = constants::alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp)*factorxp/(k1*k1)*H_co
        *(StF(k,TA,QsA,mv, BK, xA)*factorxA*StF(pminuskminusk1,TA,QsA,mv, BK, xA)*factorxA)
        *Rscale*Rscale
        *p*pscale*2.*constants::PI
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // dRxdRy
  // d2p
  // d2k
  // d2k1
    }
  }
  return 0;
}

// Integrand for the combined J/Psi integral
static int JPsiIntegrandAllNob(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobqqqtilde xx[0]
#define nobqqphi xx[1]
#define nobqqM xx[2]
#define nobqqPT xx[3]
#define nobqq4k xx[4]
#define nobqq4phik xx[5]
#define nobqq4k1 xx[6]
#define nobqq4phik1 xx[7]

  double kscale = 30.;
  double pscale = 30.;
  // Qs will be made rapidity dependent
  double BK = static_cast<params*>(userdata)->BK;
  
  double Y = static_cast<params*>(userdata)->Y;
  // double Qsp = static_cast<params*>(userdata)->Qsp;
  // double QsA = static_cast<params*>(userdata)->QsA;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + nobqqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = nobqqqtilde*qtildescale;
  double PT = nobqqPT*pscale; 
  double k = nobqq4k*kscale;
  double phik = nobqq4phik*2.*constants::PI;
  double k1 = nobqq4k1*kscale;
  double phik1 = nobqq4phik1*2.*constants::PI;

  kinPair in;
  in.M = M;
  in.PT  = PT;
  in.phi = nobqqphi*2.*constants::PI; // not the PT phi
  in.qtilde = qtilde;
  in.Y = Y;
  in.m = m;

  kinqq out = convert(in);
  
  double p = out.p;
  double phip = out.phip;
  double q = out.q;
  double phiq = out.phiq;
  double yp = out.yp;
  double yq = out.yq;
  
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double qx = q*cos(phiq); 
  double qy = q*sin(phiq);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pplusqminusk1minuskx = px+qx-kx-k1x;
  double pplusqminusk1minusky = py+qy-ky-k1y;
  double pplusqminusk1minusk = sqrt(pplusqminusk1minuskx*pplusqminusk1minuskx
                                    +pplusqminusk1minusky*pplusqminusk1minusky);
  double phi_pplusqminusk1minusk = atan2(pplusqminusk1minusky,pplusqminusk1minuskx);

  double pplusqminusk1x = px+qx-k1x;
  double pplusqminusk1y = py+qy-k1y;
  double pplusqminusk1 = sqrt(pplusqminusk1x*pplusqminusk1x+pplusqminusk1y*pplusqminusk1y);
  double phi_pplusqminusk1 = atan2(pplusqminusk1y,pplusqminusk1x);

  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/constants::roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/constants::roots;

  //!!
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else{
    double betax = PT/sqrt(M*M+PT*PT);
    double gammax = 1./sqrt(1.-betax*betax);
    // get Jacobian  
    double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));
    
    double myTA = 1.*constants::bindep_A;
    double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);

    double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

    f = constants::alphas*double(constants::Nc)*double(constants::Nc)
      /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, 0., Qsp, sizeFactor, mv, BK,xp)*factorxp/(k1*k1)*H*J
      *(StF(pplusqminusk1minusk,myTA,QsA,mv, BK, xA)*factorxA*StF(k,myTA,QsA,mv, BK, xA)*factorxA)
      *2.*constants::PI*constants::Bp
      *2.*constants::PI*constants::Bt
      *PT*pscale*2.*constants::PI
      *(2.*constants::mD-constants::mJPsi)*2.*M*(M/constants::mJPsi)*(M/constants::mJPsi)
      *qtildescale
      *2.*constants::PI
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *2.; // for p and q direction
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
  // d2R
  // d2b
  // d2PT
  // dM
  // dqtilde
  // dphi
  // d2k
  // d2k1
  
  }
  return 0;
}  

// Integrand for the combined J/Psi integral
static int JPsiIntegrandAll(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define qqR xx[0]
#define qqphiR xx[1]
#define qqb xx[2]
#define qqphib xx[3]
#define qqqtilde xx[4]
#define qqphi xx[5]
#define qqM xx[6]
#define qqPT xx[7]
#define qq4k xx[8]
#define qq4phik xx[9]
#define qq4k1 xx[10]
#define qq4phik1 xx[11]

  // double kscale = 10.;
  // double pscale = 10.;
  // double Rscale = 20./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  // double bscale = 20./constants::hbarc; 
  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; 
  // Qs will be made rapidity dependent
  double BK = static_cast<params*>(userdata)->BK;
  
  double Y = static_cast<params*>(userdata)->Y;
  // double Qsp = static_cast<params*>(userdata)->Qsp;
  // double QsA = static_cast<params*>(userdata)->QsA;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + qqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = qqqtilde*qtildescale;
  double PT = qqPT*pscale; 
  //double phiPT = qqphiPT*2.*constants::PI;
  double R = qqR*Rscale;
  double b = qqb*bscale;
  double k = qq4k*kscale;
  double phik = qq4phik*2.*constants::PI;
  double k1 = qq4k1*kscale;
  double phik1 = qq4phik1*2.*constants::PI;
  double phiR = qqphiR*2*constants::PI;
  double phib = qqphib*2*constants::PI;

  kinPair in;
  in.M = M;
  in.PT  = PT;
  in.phi = qqphi*2.*constants::PI; // not the PT phi
  in.qtilde = qtilde;
  in.Y = Y;
  in.m = m;

  kinqq out = convert(in);
  
  double p = out.p;
  double phip = out.phip;
  double q = out.q;
  double phiq = out.phiq;
  double yp = out.yp;
  double yq = out.yq;
  
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double qx = q*cos(phiq); 
  double qy = q*sin(phiq);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pplusqminusk1minuskx = px+qx-kx-k1x;
  double pplusqminusk1minusky = py+qy-ky-k1y;
  double pplusqminusk1minusk = sqrt(pplusqminusk1minuskx*pplusqminusk1minuskx
                                    +pplusqminusk1minusky*pplusqminusk1minusky);
  double phi_pplusqminusk1minusk = atan2(pplusqminusk1minusky,pplusqminusk1minuskx);

  double pplusqminusk1x = px+qx-k1x;
  double pplusqminusk1y = py+qy-k1y;
  double pplusqminusk1 = sqrt(pplusqminusk1x*pplusqminusk1x+pplusqminusk1y*pplusqminusk1y);
  double phi_pplusqminusk1 = atan2(pplusqminusk1y,pplusqminusk1x);
  double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
 
  // double Rx = R*cos(phiR);
  // double Ry = R*sin(phiR);
  // double bx = b*cos(phib);
  // double by = b*sin(phib);
  // double Rminusbx = Rx-bx;
  // double Rminusby = Ry-by;
  // double Rminusb = sqrt(Rminusbx*Rminusbx+Rminusby*Rminusby);
  // double phi_Rminusb = atan2(Rminusby,Rminusbx);
  
  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/constants::roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/constants::roots;

  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else{
    double betax = PT/sqrt(M*M+PT*PT);
    double gammax = 1./sqrt(1.-betax*betax);
    // get Jacobian  
    double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));
    
    double myTA = returnTA(Rminusb,TAclass)*constants::bdep_A; //(2.37=0.4*(208)^(1/3));
    double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);

    double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

    f = constants::alphas*double(constants::Nc)*double(constants::Nc)
      /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor, mv, BK,xp)*factorxp/(k1*k1)*H*J
      *(StF(pplusqminusk1minusk,myTA,QsA,mv, BK, xA)*factorxA*StF(k,myTA,QsA,mv, BK, xA)*factorxA)
      *R*Rscale*2.*constants::PI
      *b*bscale*2.*constants::PI
      *PT*pscale*2.*constants::PI
      *(2.*constants::mD-constants::mJPsi)*2.*M*(M/constants::mJPsi)*(M/constants::mJPsi)
      *qtildescale
      *2.*constants::PI
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *2.; // for p and q direction
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
  // d2R
  // d2b
  // d2PT
  // dM
  // dqtilde
  // dphi
  // d2k
  // d2k1
  
  }
  return 0;
}  

// Integrand for the combined J/Psi integral (no b integral) for fluctuations
static int JPsiIntegrandAllFluc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define fqqRx xx[0]
#define fqqRy xx[1]
#define fqqqtilde xx[2]
#define fqqphi xx[3]
#define fqqM xx[4]
#define fqqPT xx[5]
#define fqq4k xx[6]
#define fqq4phik xx[7]
#define fqq4k1 xx[8]
#define fqq4phik1 xx[9]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  // double bscale = 24./constants::hbarc; 
  // double Qsp = static_cast<params*>(userdata)->Qsp;
  // double QsA = static_cast<params*>(userdata)->QsA;
  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
 
  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + fqqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = fqqqtilde*qtildescale;
  double PT = fqqPT*pscale; 
  //double phiPT = qqphiPT*2.*constants::PI;

  double Rx = fqqRx*Rscale-Rscale/2.;
  double Ry = fqqRy*Rscale-Rscale/2.;
  //double bx = fqqbx*bscale-bscale/2.;
  //double by = fqqby*bscale-bscale/2.;
  //double R = fqqR*Rscale;
  //double b = fqqb*bscale;
  //double phiR = fqqphiR*2*constants::PI;
  //double phib = fqqphib*2*constants::PI;
  double k = fqq4k*kscale;
  double phik = fqq4phik*2.*constants::PI;
  double k1 = fqq4k1*kscale;
  double phik1 = fqq4phik1*2.*constants::PI;

  kinPair in;
  in.M = M;
  in.PT  = PT;
  in.phi = fqqphi*2.*constants::PI; // not the PT phi
  in.qtilde = qtilde;
  in.Y = Y;
  in.m = m;

  kinqq out = convert(in);
  
  double p = out.p;
  double phip = out.phip;
  double q = out.q;
  double phiq = out.phiq;
  double yp = out.yp;
  double yq = out.yq;

  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/constants::roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/constants::roots;

  double factorxp = pow(1.-xp,4.);
  if (xp>1.){
    factorxp = 0.;
  }
  double factorxA = pow(1.-xA,4.);
  if (xA>1.){
    factorxA = 0.;
  }

  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double qx = q*cos(phiq); 
  double qy = q*sin(phiq);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pplusqminusk1minuskx = px+qx-kx-k1x;
  double pplusqminusk1minusky = py+qy-ky-k1y;
  double pplusqminusk1minusk = sqrt(pplusqminusk1minuskx*pplusqminusk1minuskx
                                    +pplusqminusk1minusky*pplusqminusk1minusky);
  double phi_pplusqminusk1minusk = atan2(pplusqminusk1minusky,pplusqminusk1minuskx);

  double pplusqminusk1x = px+qx-k1x;
  double pplusqminusk1y = py+qy-k1y;
  double pplusqminusk1 = sqrt(pplusqminusk1x*pplusqminusk1x+pplusqminusk1y*pplusqminusk1y);
  double phi_pplusqminusk1 = atan2(pplusqminusk1y,pplusqminusk1x);
  //double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
 
  // double Rx = R*cos(phiR);
  // double Ry = R*sin(phiR);
  // double bx = b*cos(phib);
  // double by = b*sin(phib);
  // double Rminusbx = Rx-bx;
  // double Rminusby = Ry-by;
  // double Rminusb = sqrt(Rminusbx*Rminusbx+Rminusby*Rminusby);
  // double phi_Rminusb = atan2(Rminusby,Rminusbx);
  
  double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);

  // get Jacobian
  double betax = PT/sqrt(M*M+PT*PT);
  double gammax = 1./sqrt(1.-betax*betax);
  
  double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));

  //double R = sqrt(Rx*Rx+Ry*Ry);
  //double b = sqrt(bx*bx+by*by);
  //double TA = returnTA(sqrt((Rx-bx)*(Rx-bx)+(Ry-by)*(Ry-by)),TAclass);
  
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass)*constants::bdep_fluc_A;
  double Tp = returnTp2D(Rx,Ry,glauberClass)*constants::bdep_fluc_p;

  // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
  f = constants::alphas*double(constants::Nc)*double(constants::Nc)
    /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp)*factorxp/(k1*k1)*H*J
    *(StF(pplusqminusk1minusk,TA,QsA,mv, BK, xA)*factorxA*StF(k,TA,QsA,mv,BK,xA))*factorxA
    *Rscale*Rscale
    //  *bscale*bscale
    *PT*pscale*2.*constants::PI
    *(2.*constants::mD-constants::mJPsi)*2.*M*(M/constants::mJPsi)*(M/constants::mJPsi)
    *qtildescale
    *2.*constants::PI
    *k*kscale*2.*constants::PI
    *k1*kscale*2.*constants::PI
    *2.; // for p and q direction
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // dRxdRy
  //// dbxdby
  // d2PT
  // dM
  // dqtilde
  // dphi
  // d2k
  // d2k1

  return 0;
}  

// Integrand for the combined ccbar integral (no b integral) for fluctuations
static int ccBarIntegrandAllFluc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define fccRx xx[0]
#define fccRy xx[1]
#define fccq xx[2]
#define fccqphi xx[3]
#define fccp xx[4]
#define fccpphi xx[5]
#define fcc4k xx[6]
#define fcc4phik xx[7]
#define fcc4k1 xx[8]
#define fcc4phik1 xx[9]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
 
  double m = constants::mc;//static_cast<params*>(userdata)->m;
  double BK = static_cast<params*>(userdata)->BK;
 
  // scale the integration variables 
  double q = fccq*pscale; // cbar transverse momentum
  double phiq = fccqphi*2.*constants::PI; // and its angle
  double p = fccp*pscale; // c transverse momentum
  double phip = fccpphi*2.*constants::PI; // and its angle
  double k = fcc4k*kscale; // exchaged k momentum
  double phik = fcc4phik*2.*constants::PI;// and its angle
  double k1 = fcc4k1*kscale; // exchaged k1 momentum
  double phik1 = fcc4phik1*2.*constants::PI;
  double Rx = fccRx*Rscale-Rscale/2.;
  double Ry = fccRy*Rscale-Rscale/2.;

  double yq = static_cast<params*>(userdata)->yq; // c bar rapidity
  double yp = static_cast<params*>(userdata)->yp; // c rapidity

  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/constants::roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/constants::roots;

  double Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
  
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double qx = q*cos(phiq); 
  double qy = q*sin(phiq);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pplusqminusk1minuskx = px+qx-kx-k1x;
  double pplusqminusk1minusky = py+qy-ky-k1y;
  double pplusqminusk1minusk = sqrt(pplusqminusk1minuskx*pplusqminusk1minuskx
                                    +pplusqminusk1minusky*pplusqminusk1minusky);
  double phi_pplusqminusk1minusk = atan2(pplusqminusk1minusky,pplusqminusk1minuskx);

  double pplusqminusk1x = px+qx-k1x;
  double pplusqminusk1y = py+qy-k1y;
  double pplusqminusk1 = sqrt(pplusqminusk1x*pplusqminusk1x+pplusqminusk1y*pplusqminusk1y);
  double phi_pplusqminusk1 = atan2(pplusqminusk1y,pplusqminusk1x);
 
  double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);
  // double H = Hard::qqqq(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m)
  //  +Hard::qqg(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m)
  //  +Hard::gg(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);
 
  double R = sqrt(Rx*Rx+Ry*Ry);
  
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass);
  double Tp = returnTp2D(Rx,Ry,glauberClass);

  // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
  f = constants::alphas*double(constants::Nc)*double(constants::Nc)
    /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp)/(k1*k1)*H
    *(StF(pplusqminusk1minusk,TA,QsA,mv, BK, xA)*StF(k,TA,QsA,mv, BK, xA))
    *Rscale*Rscale // dRx dRy
    *k*kscale*2.*constants::PI // d2k
    *k1*kscale*2.*constants::PI // d2k1
    *p*pscale*2.*constants::PI // d2p
    *q*pscale*2.*constants::PI // d2q
    *2.; // for p and q direction
 

  return 0;
}  
// Integrand for the combined J/Psi integral
static int JPsiIntegrandNoPT(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nqqR xx[0]
#define nqqphiR xx[1]
#define nqqb xx[2]
#define nqqphib xx[3]
#define nqqqtilde xx[4]
#define nqqphi xx[5]
#define nqqM xx[6]
#define nqq4k xx[7]
#define nqq4phik xx[8]
#define nqq4k1 xx[9]
#define nqq4phik1 xx[10]

  double kscale = 30.;
  double Rscale = 2./constants::hbarc; // choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; 
  // Qs will be made rapidity dependent
  double Qsp = static_cast<params*>(userdata)->Qsp;
  double QsA = static_cast<params*>(userdata)->QsA;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  double m = constants::mc;//static_cast<params*>(userdata)->m;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  TAInt* TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  double BK = static_cast<params*>(userdata)->BK;
 

  // scale the integration variables 
  double M = constants::mJPsi + nqqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = nqqqtilde*qtildescale;
  double R = nqqR*Rscale;
  double b = nqqb*bscale;
  double k = nqq4k*kscale;
  double phik = nqq4phik*2.*constants::PI;
  double k1 = nqq4k1*kscale;
  double phik1 = nqq4phik1*2.*constants::PI;
  double phiR = nqqphiR*2*constants::PI;
  double phib = nqqphib*2*constants::PI;


  kinPair in;
  in.M = M;
  in.PT  = PT;
  in.phi = nqqphi*2.*constants::PI; // not the PT phi
  in.qtilde = qtilde;
  in.Y = Y;
  in.m = m;

  kinqq out = convert(in);
  
  double p = out.p;
  double phip = out.phip;
  double q = out.q;
  double phiq = out.phiq;
  double yp = out.yp;
  double yq = out.yq;
  
  // get sums of vectors
  double px = p*cos(phip); 
  double py = p*sin(phip);
  double qx = q*cos(phiq); 
  double qy = q*sin(phiq);
  double kx = k*cos(phik);
  double ky = k*sin(phik);
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
  double pplusqminusk1minuskx = px+qx-kx-k1x;
  double pplusqminusk1minusky = py+qy-ky-k1y;
  double pplusqminusk1minusk = sqrt(pplusqminusk1minuskx*pplusqminusk1minuskx
                                    +pplusqminusk1minusky*pplusqminusk1minusky);
  double phi_pplusqminusk1minusk = atan2(pplusqminusk1minusky,pplusqminusk1minuskx);

  double pplusqminusk1x = px+qx-k1x;
  double pplusqminusk1y = py+qy-k1y;
  double pplusqminusk1 = sqrt(pplusqminusk1x*pplusqminusk1x+pplusqminusk1y*pplusqminusk1y);
  double phi_pplusqminusk1 = atan2(pplusqminusk1y,pplusqminusk1x);
  double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));

  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/constants::roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/constants::roots;

 
  // double Rx = R*cos(phiR);
  // double Ry = R*sin(phiR);
  // double bx = b*cos(phib);
  // double by = b*sin(phib);
  // double Rminusbx = Rx-bx;
  // double Rminusby = Ry-by;
  // double Rminusb = sqrt(Rminusbx*Rminusbx+Rminusby*Rminusby);
  // double phi_Rminusb = atan2(Rminusby,Rminusbx);
  
  double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);

  // get Jacobian
  double betax = PT/sqrt(M*M+PT*PT);
  double gammax = 1./sqrt(1.-betax*betax);
  double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));

  double myTA = returnTA(Rminusb,TAclass)*constants::bdep_A; //TA(Rminusb);

  f = constants::alphas*double(constants::Nc)*double(constants::Nc)
    /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, R, Qsp, sizeFactor, mv, BK,xp)/(k1*k1)*H*J
    *(StF(pplusqminusk1minusk,myTA,QsA,mv, BK, xA)*StF(k,myTA,QsA,mv, BK, xA))
    *R*Rscale*2.*constants::PI
    *b*bscale*2.*constants::PI
    *(2.*constants::mD-constants::mJPsi)*2.*M*(M/constants::mJPsi)*(M/constants::mJPsi)
    *qtildescale
    *2.*constants::PI
    *k*kscale*2.*constants::PI
    *k1*kscale*2.*constants::PI
    *2.; // for p and q direction
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // dM
  // dqtilde
  // dphi
  // d2k
  // d2k1
  
  //don't forget F, it is one at the moment.

  return 0;
}  

static int GluonsNoB(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobk xx[0]
#define nobphik xx[1]
#define nobp xx[2]
#define nobphi xx[3]

  double kscale = 30.;
  double pscale = 30.;
  
  double Qsp;// = static_cast<params*>(userdata)->Qsp;
  double QsA;// = static_cast<params*>(userdata)->QsA;
  double lambda = static_cast<params*>(userdata)->lambda;
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
 
  double xp = (nobp*pscale+lambda)*exp(Y)/constants::roots;
  double xA = (nobp*pscale+lambda)*exp(-Y)/constants::roots;
  Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

  double TA = 1.*constants::bdep_A; // To avoid impact parameter dependence. We also set R=0 inside Phip for the same purpose
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    f = constants::alphas/constants::CF/(nobp*pscale+lambda)/(nobp*pscale+lambda)/pow((2*constants::PI*constants::PI),3.)
      *Phip(nobk*kscale, 0, Qsp, sizeFactor, mv, BK, xp)*Phit(sqrt((nobp*pscale+lambda)*(nobp*pscale+lambda) + nobk*nobk*kscale*kscale - 2.*(nobp*pscale+lambda)*nobk*kscale*cos((nobphi-nobphik)*2.*constants::PI)), TA, QsA, mv, BK, xA)
      *2.*constants::PI*nobk*kscale*kscale  //kdkdphik
      *2.*constants::PI*constants::Bt  //R-integral
      *2.*constants::PI*constants::Bp  // b-integral
      *2.*constants::PI*pscale*(nobp*pscale+lambda); //pdpdphip
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
    return 1;
}

// Integrand for integral over everything but |p|
static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define gk xx[0]
#define gphik xx[1]
#define gR xx[2]
#define gphiR xx[3]
#define gb xx[4]
#define gphib xx[5]
#define gphi xx[6]

  double kscale = 15.;
  double Rscale = 10./constants::hbarc;
  double p = static_cast<params*>(userdata)->pe;
  double Qsp = static_cast<params*>(userdata)->Qsp;
  double QsA = static_cast<params*>(userdata)->QsA;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double BK = static_cast<params*>(userdata)->BK;
 
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  double TA = returnTA(sqrt(max(gR*Rscale*gR*Rscale + gb*gb*Rscale*Rscale - 2.*gR*gb*Rscale*Rscale*cos((gphiR - gphib)*2.*constants::PI),0.)),TAclass);

  // determine x values

  f = constants::alphas/constants::CF/(p)/(p)/pow((2*constants::PI*constants::PI),3.)
    *Phip(gk*kscale, gR*Rscale, Qsp, sizeFactor, mv, BK,0)*Phit(sqrt(p*p + gk*gk*kscale*kscale - 2.*p*gk*kscale*cos((gphi - gphik)*2.*constants::PI)), TA, QsA, mv, BK, 0)
    *2.*constants::PI*gk*kscale*kscale  //kdkdphik
    *2.*constants::PI*gR*Rscale*Rscale  //RdRdphiR
    *2.*constants::PI*gb*Rscale*Rscale;  //bdbdphib
  return 0;
}

static int FullIntegrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define gp xx[7]

  //cout << "xx[0] = " << xx[0] << endl;
  if (xx[0]>1.){
    cout << "xx[0] = " << xx[0] << endl;
    cout << "xx[7] = " << xx[7] << endl;
    cout << "f=" << f << endl;
 }

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 2./constants::hbarc;
  double bscale = 12./constants::hbarc; //!!
  double Qsp;// = static_cast<params*>(userdata)->Qsp;
  double QsA;// = static_cast<params*>(userdata)->QsA;
  double lambda = static_cast<params*>(userdata)->lambda;
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
 
  double xp = (gp*pscale+lambda)*exp(Y)/constants::roots;
  double xA = (gp*pscale+lambda)*exp(-Y)/constants::roots;
  double factorxA = pow(1.-xA,4.);
  double factorxp = pow(1.-xp,4.);
  if (xp>1.){
    f = 0.;
  }
  else if (xA>1.){
    f = 0.;
  }
  else{
    
    Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    double TA = returnTA(sqrt(max(gR*Rscale*gR*Rscale + gb*gb*bscale*bscale - 2.*gR*gb*Rscale*bscale*cos((gphiR - gphib)*2.*constants::PI),0.)),TAclass)*constants::bdep_A; //(2.37=0.4*(208)^(1/3));
    
    f = constants::alphas/constants::CF/(gp*pscale+lambda)/(gp*pscale+lambda)/pow((2*constants::PI*constants::PI),3.)
      *Phip(gk*kscale, gR*Rscale, Qsp, sizeFactor, mv, BK,xp)*Phit(sqrt((gp*pscale+lambda)*(gp*pscale+lambda) + gk*gk*kscale*kscale - 2.*(gp*pscale+lambda)*gk*kscale*cos((gphi - gphik)*2.*constants::PI)), TA, QsA, mv, BK, xA)
      *2.*constants::PI*gk*kscale*kscale  //kdkdphik
      *2.*constants::PI*gR*Rscale*Rscale  //RdRdphiR
      *2.*constants::PI*gb*bscale*bscale  //bdbdphib
      *2.*constants::PI*pscale*(gp*pscale+lambda); //pdpdphip
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 1;
}


static int FullIntegrandFluc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define fgk xx[0]
#define fgphik xx[1]
#define fgRx xx[2]
#define fgRy xx[3]
#define fgphi xx[4]
#define fgp xx[5]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc;
  //double bscale = 24./constants::hbarc;
  //double Rscale = 1./constants::hbarc;
  //double bscale = 4./constants::hbarc;
  double Rx = fgRx*Rscale-Rscale/2.;
  double Ry = fgRy*Rscale-Rscale/2.;
  //double bx = fgbx*bscale-bscale/2.;
  //double by = fgby*bscale-bscale/2.;

  double BK = static_cast<params*>(userdata)->BK;

  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  //double Qsp = static_cast<params*>(userdata)->Qsp;
  //double QsA = static_cast<params*>(userdata)->QsA;
  double lambda = static_cast<params*>(userdata)->lambda;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double Y = static_cast<params*>(userdata)->Y;
  double Qsp=0.;
  double QsA=0.;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double R = sqrt(Rx*Rx+Ry*Ry);
  double b = sqrt(bx*bx+by*by);
  //double TA = returnTA(sqrt((Rx-bx)*(Rx-bx)+(Ry-by)*(Ry-by)),TAclass);
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass)*constants::bdep_fluc_A;
  double Tp = returnTp2D(Rx,Ry,glauberClass)*constants::bdep_fluc_p;

  double xp = (fgp*pscale+lambda)*exp(Y)/constants::roots;
  double xA = (fgp*pscale+lambda)*exp(-Y)/constants::roots;
  
  double factorxA = pow(1.-xA,4.);
  double factorxp = pow(1.-xp,4.);
  if (xp>1.){
    f = 0.;
  }
  else if (xA>1.){
    f = 0.;
  }
  else{
    
    Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
    f = constants::alphas/constants::CF/(fgp*pscale+lambda)/(fgp*pscale+lambda)/pow((2*constants::PI*constants::PI),3.)
      *PhipFluc(fgk*kscale, Tp, Qsp, sizeFactor, mv, BK, xp)*factorxp*Phit(sqrt((fgp*pscale+lambda)*(fgp*pscale+lambda) + fgk*fgk*kscale*kscale - 2.*(fgp*pscale+lambda)*fgk*kscale*cos((fgphi - fgphik)*2.*constants::PI)), TA, QsA, mv, BK, xA)*factorxA
      *2.*constants::PI*fgk*kscale*kscale  //kdkdphik
      *Rscale*Rscale  //dRxdRy
      //    *bscale*bscale  //bdxdby
      *2.*constants::PI*pscale*(fgp*pscale+lambda); //pdpdphip
  }
  
  return 1;
}


// Integrand for the full 9D integral
static int HadronIntegrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define z xx[6]

  double Kch = 1.; //irrelevant for the ratio
  double sigmainel = 1.; //irrelevant for the ratio
  double mh = 0.3; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc;
  double bscale = 24./constants::hbarc;
  //double Rscale = 1./constants::hbarc;
  //double bscale = 4./constants::hbarc;
  double Rx = fgRx*Rscale-Rscale/2.;
  double Ry = fgRy*Rscale-Rscale/2.;
  //double bx = fgbx*bscale-bscale/2.;
  //double by = fgby*bscale-bscale/2.;

  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  //double Qsp = static_cast<params*>(userdata)->Qsp;
  //double QsA = static_cast<params*>(userdata)->QsA;
  double lambda = static_cast<params*>(userdata)->lambda;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double Y = static_cast<params*>(userdata)->Y;
  double eta = Y;
  double Qsp=0.;
  double QsA=0.;

  double BK = static_cast<params*>(userdata)->BK;
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double R = sqrt(Rx*Rx+Ry*Ry);
  double b = sqrt(bx*bx+by*by);
  //double TA = returnTA(sqrt((Rx-bx)*(Rx-bx)+(Ry-by)*(Ry-by)),TAclass);
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass)*constants::bdep_fluc_A;
  double Tp = returnTp2D(Rx,Ry,glauberClass)*constants::bdep_fluc_p;
  
  double p = (fgp*pscale+lambda);
  
  double pg = p/z;
  
  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  
  double yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))/((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/constants::roots;
  double xA = pg*exp(-yg)/constants::roots;
  
  if (xp>1.){
    f = 0.;
  }
  else if (xA>1.){
    f = 0.;
  }
  else if (pg>30.){
    f = 0.;
  }    
  else{
    double factorxA = pow(1.-xA,4.);
    double factorxp = pow(1.-xp,4.);
    
    Qsp = constants::prefactor*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    QsA = constants::prefactor*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
    f = Dh/z/z*J* 
      constants::alphas/constants::CF/(pg)/(pg)/pow((2*constants::PI*constants::PI),3.)
      *PhipFluc(fgk*kscale, Tp, Qsp, sizeFactor, mv, BK, xp)*factorxp*Phit(sqrt(pg*pg + fgk*fgk*kscale*kscale - 2.*pg*fgk*kscale*cos((fgphi - fgphik)*2.*constants::PI)), TA, QsA, mv, BK, xA)*factorxA
      *2.*constants::PI*fgk*kscale*kscale  //kdkdphik
      *Rscale*Rscale  //dRxdRy
      //    *bscale*bscale  //bdxdby
      *2.*constants::PI*pscale*p; //pdpdphip
  }
  
  return 1;
}

double Qsp(double pT, double roots, double y){
  //return pow((0.0003*roots/pT/exp(-y)),0.288/2.);
  return pow((constants::x0*roots/pT/exp(-y)),constants::lambdaSpeedp/2.);
  //return pow((0.00008*roots/pT/exp(-y)),0.9/2.);
  //return 0.5;
}

double QsA(double pT, double roots, double y){
  //  return sqrt(0.4*pow(208.,(1./3.))*pow(Qsp(pT,roots,y),2.));
  return pow((constants::x0*roots/pT/exp(-y)),constants::lambdaSpeedA/2.);
}

// Main program
int main(int argc, char *argv[]) {
  // MPI things
  int rank=0;
  int size=1;
  
  // Options
  int readTable = 0;
  int useFluc = 0;
  int Nevents = 1;
  int NRQCD = 0;
  int BK = 0;

  std::vector <std::string> sources;
  std::string destination;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--readTable") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        readTable = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--readTable option requires one argument, 0 or 1." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--fluctuations") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        useFluc = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--fluctuations option requires one argument, 0 or 1." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Nevents") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        Nevents = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Nevents option requires one argument, an integer >=1." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--NRQCD") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        NRQCD = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--NRQCD option requires one argument, 0 for ICEM or 1 for NRQCD." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--BK") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        BK = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--BK option requires one argument, 0 for MV with Qs(x) or 1 for (at the moment fake) BK." << std::endl;
        return 1;
      }  
    }
  }
    
  // initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
  MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes

 
  int h5Flag = 0;
  pretty_ostream messenger;

 
  if (rank==0){
    display_logo();
    cout << "Options: read MV dipole from file yes(1)/no(0)= " << readTable << ", fluctuations on(1)/off(0) = " << useFluc << ", Number of events = " << Nevents << endl;
    messenger.flush("info");
  }
  
  long int seed = time(NULL)+rank*100000;
  //long int seed = 1;
 
  Parameters *Glauber_param;
  Glauber_param = new Parameters();
  Glauber_param->setParameters();

  Random *random;
  random = new Random();
  random->init_genrand64(seed);

  Glauber *glauber;
  glauber = new Glauber(Glauber_param);
  glauber->init(random);
  glauber->makeNuclei(random, constants::Bp);

  TAInt *TAclass;
  TAclass = new TAInt();
  TAclass->computeTAIntegral();

  // Make the MV table 
  MV *mv;
  mv = new MV(BK);

  if (readTable == 0){
    if(BK==0){
      mv->computePhip();
      mv->writeTable();
      //mv->writeTableText();
    }
    else{
      mv->computePhipBK();
      mv->writeTableBK();
    }
  }
  else if (readTable == 1){
    if(BK==0){
      mv->readTable();
    }
    else{
      mv->readTableBK();
    }
  }
  else{
    cerr << "Unknown option readTable = " << readTable << ". Options are 0 or 1." << endl;
    exit(0);
  }

  //  messenger.flush("info");

  //  cout << "Phip=" << mv->Phip(0.1,1.,1.) << endl;

  ////double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);
 
  // Cuba's parameters for integration
  int NDIM = 9;
  int NCOMP = 1;
  const long long int NVEC = 1;
  //  double EPSREL = 5e-4;
  //  double EPSREL = 5e-4;
  //double EPSABS = 1e-15;
  double EPSREL = 5e-4;
  double EPSABS = 1e-15;
  int VERBOSE = 0;
  int LAST = 4;
  const long long int MINEVAL = 0;
  
  /// put the large number back in !!! 
  //const long long int MAXEVAL = 5000000000;
  const long long int MAXEVAL =   50000000;
  int KEY = 0;
  
  //vegas
  int SEED = time(NULL)+rank*100000;
  const long long int NSTART = 10000000;
  const long long int NINCREASE = 1000000;
  const long long int NBATCH = 1000;
  int GRIDNO = 0;

  //suave
  int NNEW = 2000;
  int NMIN = 2;
  int FLATNESS = 25;

  //divonne
  int KEY1 = 7;
  int KEY2 = 7;
  int KEY3 = 1;
  int MAXPASS = 1000;
  double BORDER = 0.;
  double MAXCHISQ = 10.;
  double MINDEVIATION = 0.25;
  int NGIVEN = 0;
  int LDXGIVEN = 0;
  int NEXTRA = 0;
 
  int comp, nregions, fail;
  long long int neval;

  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  params data;

  //  double QspPre = 0.43; // prefactors for scaling
  //  double QsAPre = 0.43; // prefactors for scaling
  double QspPre = constants::prefactor; // prefactors for scaling
  double QsAPre = constants::prefactor; // prefactors for scaling

  double inQsp;
  double inQsA;

  double Y_g = 0.;
  double Y_fwd = 3.;
  double Y_bck= -3.8; //use minus sign
 
  data.Y = Y_g;

  data.BK = BK;

  if(useFluc ==0){
    inQsp = QspPre*Qsp(0.8,constants::roots,-Y_g);
    inQsA = QsAPre*QsA(0.8,constants::roots,Y_g);
  }
  else{
    //    inQsp = QspPre*Qsp(0.8,constants::roots,0.);
    //inQsA = QsAPre*Qsp(0.8,constants::roots,0.);
    inQsp = QspPre*Qsp(0.8,constants::roots,-Y_g);
    inQsA = QsAPre*QsA(0.8,constants::roots,Y_g);
  }

  data.PT = 0.; // dummy for now
  data.pe = 0.; // dummy for now
  data.k = 0.;  // dummy for now
  data.m = 0.;
  data.Qs = 0.; // Saturation scale in GeV - not used anymore
  data.lambda = 0.05; // Infrared cutoff on p integral in GeV (50 MeV according to https://arxiv.org/pdf/1812.01312.pdf)
  data.Qsp = inQsp; // midrapidity proton Saturation scale in GeV
  data.QsA = inQsA; // midrapidity Pb Saturation scale in GeV
  data.mv = mv; // MV class
  data.TAclass = TAclass; // TA class
  data.glauberClass = glauber; // Glauber class
  data.protonSizeFactor = 1.; // allows to make proton larger at small x and smaller at large x (not used, leave at one)
 
  if (rank==0){
  cout << "Qsp(y=0) = " << inQsp << endl;
  cout << "QsA(y=0) = " << inQsA << endl;
  }

  double inQsp_fwd;
  double inQsA_fwd;

  if(useFluc == 0){
    inQsp_fwd = QspPre*Qsp(3,constants::roots,-Y_fwd);
    inQsA_fwd = QsAPre*QsA(3,constants::roots,Y_fwd);
  }
  else{
    //    inQsp_fwd = QspPre*Qsp(3,constants::roots,-3.);
    //inQsA_fwd = QsAPre*Qsp(3,constants::roots,3.);
    inQsp_fwd = QspPre*Qsp(3.,constants::roots,-Y_fwd);
    inQsA_fwd = QsAPre*QsA(3.,constants::roots,Y_fwd);
  }

  if (rank==0){
    cout << "Qsp(y=" << Y_fwd << ") = " << inQsp_fwd << endl;
    cout << "QsA(y="<< Y_fwd << ") = " << inQsA_fwd << endl;
  }

  double inQsp_bck;
  double inQsA_bck;

  if(useFluc == 0){
    inQsp_bck = QspPre*Qsp(2.7,constants::roots,-Y_bck);
    inQsA_bck = QsAPre*QsA(2.7,constants::roots,Y_bck);
  }
  else{
    //   inQsp_bck = QspPre*Qsp(2.7,constants::roots,3.8);
    //inQsA_bck = QsAPre*Qsp(2.7,constants::roots,-3.8);
    inQsp_bck = QspPre*Qsp(2.7,constants::roots,-Y_bck);
    inQsA_bck = QsAPre*QsA(2.7,constants::roots,Y_bck);
  }

  if (rank==0){
  cout << "Qsp(y=" << Y_bck << ") = " << inQsp_bck << endl;
  cout << "QsA(y=" << Y_bck << ") = " << inQsA_bck << endl;
  }
  // // test the interpolation routine
  // for(int i=0;i<100;i++){
  //   double myR = double(i)/10./constants::hbarc;
  //   cout << returnTA(myR,TAclass) <<  " " << TA(myR) << endl;
  // }


 
  double gresult;
  double gerror;

  double hresult;
  double herror;

  double JPsi2result;
  double JPsi2error;

  double JPsi2result_cs;
  double JPsi2error_cs;

  double JPsi2result_co;
  double JPsi2error_co;

  double JPsi2result2;
  double JPsi2error2;

  //test StF
 // for (int ik=0; ik<1000; ik++){
  //  double k = ik*0.01;
    //double StF(double k, double TA, double Qs, MV *mv, int BK, double x)
  //  cout << k << " " << StF(k, 1., 0.4949, mv, 0, 0.0001)  << endl;
 // }


  // Now compute midrapidity gluons
  if(useFluc == 0){
    cout << "For b integrated results obtained in this mode (no fluctuations) all results are cross sections, that need to be divided by the total inelastic cross section (in p+Pb) to get particle numbers." << endl;
    // Run 8D Vegas integration
    // printf("Midrapidity gluon: %.8f +- %.8f\t\n", gresult, gerror);

    //data.Qsp = inQsp_fwd; // forward proton Saturation scale in GeV
    //data.QsA = inQsA_fwd; // forward Pb Saturation scale in GeV
    //data.Y = Y_fwd;
    double ptmin = log10(0.01);
    double ptmax = log10(5);
    int npoints = 20;
    
    double step = (ptmax-ptmin)/npoints;
    // stringstream strfilename;
    // strfilename << "djpsi_gluons_bk_nrqcd_nob.dat";
    // string filename;
    // filename = strfilename.str();
    // fstream fout(filename.c_str(), ios::app);
    // for (int nip=0; nip<=npoints; nip++){ 
    // data.Y=ymin + nip*step;

     
    // NDIM = 4;
    // llVegas(NDIM, NCOMP, GluonsNoB, &data, NVEC,
    //           EPSREL, EPSABS, VERBOSE, SEED,
    //           MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    //           GRIDNO, NULL, NULL,
    //           &neval, &fail, integral, error, prob);

    // gresult = (double)integral[0];
    // gerror = (double)error[0];          
    // // JPsi cross section
    // if(NRQCD==1){

    //   NDIM = 8;
    //   llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsNob, &data, NVEC,
    //           EPSREL, EPSABS, VERBOSE, SEED,
    //           MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    //           GRIDNO, NULL, NULL,
    //           &neval, &fail, integral, error, prob);

    //   JPsi2result_cs = (double)integral[0];
    //   JPsi2error_cs = (double)error[0];

    //   NDIM = 6;
    //   llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoNob, &data, NVEC,
    //           EPSREL, EPSABS, VERBOSE, SEED,
    //           MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    //           GRIDNO, NULL, NULL,
    //           &neval, &fail, integral, error, prob);

    //   JPsi2result_co= (double)integral[0];
    //   JPsi2error_co = (double)error[0];
    // }
    // else{   
 
    //   NDIM = 8;
    //   llVegas(NDIM, NCOMP, JPsiIntegrandAllNob, &data, NVEC,
    //           EPSREL, EPSABS, VERBOSE, SEED,
    //           MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    //           GRIDNO, NULL, NULL,
    //           &neval, &fail, integral, error, prob);
    //   JPsi2result = (double)integral[0];
    //   JPsi2error = (double)error[0];
    // }
    stringstream strfilename;
    strfilename << "avpt_den_djpsi_bk_nrqcd_pt_nofluc_bind_bwd.dat";
    string filename;
    filename = strfilename.str();
    fstream fout(filename.c_str(), ios::app);
    data.Y = Y_bck;
    

   // for (int nip=0; nip<=npoints; nip++){ 
   //     data.PT= pow(10,ptmin + nip*step);

     
     //NDIM = 8;
     //llVegas(NDIM, NCOMP, FullIntegrand, &data, NVEC,
     //          EPSREL, EPSABS, VERBOSE, SEED,
      //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
      //         GRIDNO, NULL, NULL,
      //         &neval, &fail, integral, error, prob);

     //gresult = (double)integral[0];
    // gerror = (double)error[0];          
    // JPsi cross section
    if(NRQCD==1){

      NDIM = 12;
      llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsNob, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);

      JPsi2result_cs = (double)integral[0];
      JPsi2error_cs = (double)error[0];

      NDIM = 10;
      llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoNob, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);

      JPsi2result_co= (double)integral[0];
      JPsi2error_co = (double)error[0];
    }
    else{   
      NDIM = 8;
      llVegas(NDIM, NCOMP, JPsiIntegrandAll, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);
      JPsi2result = (double)integral[0];
      JPsi2error = (double)error[0];
    }


   // Print the result
   // printf("Forward JPsi: %.8e +- %.8e\t\n", JPsi2result_cs, JPsi2error);

   // data.Qsp = inQsp_bck; // forward proton Saturation scale in GeV
   // data.QsA = inQsA_bck; // forward Pb Saturation scale in GeV
   // data.Y = Y_bck;
    
    //if(NRQCD==1){
    //  NDIM = 12;
    //  llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCs, &data, NVEC,
    //          EPSREL, EPSABS, VERBOSE, SEED,
     //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
     //         GRIDNO, NULL, NULL,
     //         &neval, &fail, integral, error, prob);

  //    JPsi2result2 = (double)integral[0];
    //  JPsi2error2 = (double)error[0];

    //  NDIM = 10;
    //  llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCo, &data, NVEC,
    //          EPSREL, EPSABS, VERBOSE, SEED,
    //          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    //          GRIDNO, NULL, NULL,
    //          &neval, &fail, integral, error, prob);

    //  JPsi2result2 += (double)integral[0];
    //  JPsi2error2 += (double)error[0];
    //}
    //else{
    //  NDIM = 12;
     // llVegas(NDIM, NCOMP, JPsiIntegrandAll, &data, NVEC,
    //          EPSREL, EPSABS, VERBOSE, SEED,
    //          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
     //         GRIDNO, NULL, NULL,
     //         &neval, &fail, integral, error, prob);
     // JPsi2result2 = (double)integral[0];
    //  JPsi2error2 = (double)error[0];
   // }

   // printf("Backward JPsi: %.8e +- %.8e\t\n", JPsi2result2, JPsi2error2);

    if(NRQCD==1){
      cout << setprecision(10)  << JPsi2result_co << " " << JPsi2error_co << " " << JPsi2result_cs << " " << JPsi2error_cs << " " << data.Y << endl;
    }
    else{
      cout << setprecision(10) << gresult << " " << gerror << " " << JPsi2result << " " << JPsi2error << " "  << data.Y << endl;
    }
    //stringstream strfilename;
    //strfilename << "output_bIntegrated.dat";
   // string filename;
   // filename = strfilename.str();
   // fstream fout(filename.c_str(), ios::app);
    
    fout << std::scientific << setprecision(5) << JPsi2result_co << " " << JPsi2error_co << " " << JPsi2result_cs << " " << JPsi2error_cs << " " << data.Y << endl;
   // }
    fout.close();


  }
  else{
    for (int ni=0; ni<Nevents; ni++){
      // Run Vegas integration with fluctuations
      // Make a new target
      glauber->makeNuclei(random, constants::Bp);

      // Sample b
      double bmin = 0.;
      double bmax = 10.;
      
      double xb =
        random->genrand64_real1(); // uniformly distributed random variable
      double b = sqrt((bmax * bmax - bmin * bmin) * xb + bmin * bmin);
      double phib = 2.*constants::PI*random->genrand64_real1();
      
      data.bx = b*cos(phib);
      data.by = b*sin(phib);
      
      cout << "Using impact parmater b=" << b << " [fm], phib=" << phib << endl;

      // fluctuate the proton Qs:
      // it is included when using quarks in the proton
      double QspFac = 1.;
      //    QspFac = sqrt((exp(random->Gauss(0, 0.5))) /
      //              std::exp(0.5 * 0.5 / 2.0)); //
      
      //      cout << "QspFac=" << QspFac << endl;
      
      // Do gluons:

      data.Qsp = inQsp*QspFac; // forward proton Saturation scale in GeV
      data.QsA = inQsA; // forward Pb Saturation scale in GeV
      data.Y = Y_g;

      NDIM = 6;
      llVegas(NDIM, NCOMP, FullIntegrandFluc, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);
      
      // Print the result
      gresult = (double)integral[0];
      gerror = (double)error[0];
      printf("Midrapidity gluon (fluc): %.8f +- %.8f\t\n", gresult, gerror);
            
      if(gresult<1.){
        cout << "Gluon number < 1, skipping event" << endl;
        continue;
      }

      // do hadrons next (one more (z) integral)
      NDIM = 7;
      llVegas(NDIM, NCOMP, HadronIntegrand, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);
      
      // Print the result
      hresult = (double)integral[0];
      herror = (double)error[0];
      printf("Hadrons (fluc): %.8f +- %.8f\t\n", hresult, herror);
            
      if(gresult<1.){
        cout << "Gluon number < 1, skipping event" << endl;
        continue;
      }

      // Gluons and hadrons done, her comes J/Psi
      
      data.Qsp = inQsp_fwd*QspFac; // forward proton Saturation scale in GeV
      data.QsA = inQsA_fwd; // forward Pb Saturation scale in GeV
      data.Y = Y_fwd;

      if(NRQCD==1){
        NDIM = 10;
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        
        JPsi2result = (double)integral[0];
        JPsi2error = (double)error[0];

        NDIM = 8;
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        
        JPsi2result += (double)integral[0];
        JPsi2error += (double)error[0]; //should add in quadrature
      }
      else{
        NDIM = 10;
        llVegas(NDIM, NCOMP, JPsiIntegrandAllFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        JPsi2result = (double)integral[0];
        JPsi2error = (double)error[0];
      }
      // Print the result
      printf("Forward JPsi: %.8e +- %.8e\t\n", JPsi2result, JPsi2error);   
      
      data.Qsp = inQsp_bck*QspFac; // forward proton Saturation scale in GeV
      data.QsA = inQsA_bck; // forward Pb Saturation scale in GeV
      data.Y = Y_bck;
   
      if(NRQCD==1){
        NDIM = 10;   
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        
        JPsi2result2 = (double)integral[0];
        JPsi2error2 = (double)error[0];

        NDIM = 8;
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        
        JPsi2result2 += (double)integral[0];
        JPsi2error2 += (double)error[0]; //should add in quadrature
      }
      else{
        NDIM = 10;
        llVegas(NDIM, NCOMP, JPsiIntegrandAllFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        JPsi2result2 = (double)integral[0];
        JPsi2error2 = (double)error[0];
      }
      // Print the result
      printf("Backward JPsi: %.8e +- %.8e\t\n", JPsi2result2, JPsi2error2);   

      double TA = returnTA2D(-data.bx,-data.by,glauber);

      cout << setprecision(10) << gresult << " " << gerror << " " << JPsi2result 
           << " " << JPsi2error << " " << JPsi2result2 << " " << JPsi2error2 
           << " " << sqrt(data.bx*data.bx+data.by*data.by) 
           << " " << TA << endl;
      
      stringstream strfilename;
      strfilename << "output_g_" << rank << ".dat";
      string filename;
      filename = strfilename.str();
      fstream fout(filename.c_str(), ios::app);
      
      fout << std::scientific << setprecision(5) << gresult << " " << gerror << " " << JPsi2result 
           << " " << JPsi2error << " " << JPsi2result2 << " " << JPsi2error2 << " " 
           << sqrt(data.bx*data.bx+data.by*data.by) 
           << " " << TA << endl;
      fout.close();

      stringstream strfilenameh;
      strfilenameh << "output_h_" << rank << ".dat";
      string filenameh;
      filenameh = strfilenameh.str();
      fstream fouth(filenameh.c_str(), ios::app);
      
      fouth << std::scientific << setprecision(5) << hresult << " " << herror << " " << JPsi2result 
            << " " << JPsi2error << " " << JPsi2result2 << " " << JPsi2error2 << endl;
      fouth.close();

    }
  }
  
  cout << " - - - - - - - - - - - - - - - - - " << endl;

  // // Integrate 11D to get PT-spectrum
  // int ppoints = 30; // Points in |p| to compute
  // double pstep = 0.25; // Step width in |p|

  // NDIM = 11;
  // int runs = 1;
  // for (int r=0; r<runs; r++){
  //   for (int i=0; i<=ppoints; i++){
  //     data.PT = 0.1+i*(pstep);
  //     SEED = time(NULL)+r*10000;
      
  //     llVegas(NDIM, NCOMP, JPsiIntegrandNoPT, &data, NVEC,
  //           EPSREL, EPSABS, VERBOSE, SEED,
  //           MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
  //           GRIDNO, NULL, NULL,
  //           &neval, &fail, integral, error, prob);
      
  //     JPsi2result = (double)integral[0];
  //     JPsi2error = (double)error[0];
  //     printf("%.3f \t \t%.8e \t%.8e\n", data.PT, JPsi2result, JPsi2error);
  //   }
  // }
    
  delete Glauber_param;
  delete random;
  delete mv;
  delete TAclass;
  delete glauber;

  MPI_Finalize();

  return 1;
}


void display_logo() {
  cout << endl;
  cout << "- compute JPsi production with fluctuations -------------------------------------------------------------------" << endl;
  cout << endl;
}
