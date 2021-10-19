
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
#include <random>
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
#include "kkp.h"

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

void display_logo();

namespace constants {
  const double PI = 3.14159265358979323846;
  const double Nc = 3.;
  const double hbarc = 0.1973269804;
  const double CA = double(Nc);
  const double CF = (double(Nc)*double(Nc) - 1.)/(2.*double(Nc));
  //  const double Bp = 4.; // !!!! <==== there is a Bp defined in MV.cpp as well - make sure they are the same || use 4 for non-fluc, 3 for fluc
  //  const double Bq = 0.3; // size of hot spots
  //  const double Bt = 1061; // == R_t = 1.1 * A^{1/3} fm ~6.5 fm
  const double mD = 1.867;
  //const double mc = 1.275; //vary? 1.4?
  const double mJPsi = 3.096916;
  const double mc = mJPsi/2.; 
  const double x0 = 0.000041;
  //const double lambdaSpeed = 0.277;
  //const double x0 = 0.00005;

  const double lambdaSpeedp = 0.277;
  const double lambdaSpeedA = 0.277;
  const double prefactorp = 0.52;
  const double prefactorA = 0.52;

  //const double lambdaSpeedp = 0.17;
  //const double lambdaSpeedA = 0.17;
  //const double prefactor = 0.48;
  //const double prefactor = 0.7;

  //const double roots = 8160.;
  const double ldme_singlet = 1.16/2./Nc; // GeV^3
  //const double ldme_singlet = 1.32; // GeV^3 arXiv:1009.5662
  const double ldme_octet_s10 = 0.089; // +- 0.0098 GeV^3 //1201.2675
  const double ldme_octet_s13 = 0.0030; // +- 0.0012 GeV^3 //1201.2675.
  const double ldme_octet_p3j = 0.0056*mc*mc; // (+- 0.0021 GeV^3) [GeV^5] //1201.2675.

  //const double ldme_octet_s10 = 0.0792; // lower limit
  //const double ldme_octet_s13 = 0.0018; // lower limit
  //const double ldme_octet_p3j = 0.0035*mc*mc; // lower limit
 
 //const double ldme_octet_s10 = 0.045; // +/- 0.0072 GeV^3 // arXiv:1009.5662
  //const double ldme_octet_s13 = 0.00312; // +- 0.00093 GeV^3 // arXiv:1009.5662 
  //const double ldme_octet_p3j = -0.0121; // (+- 0.0035 GeV^5) [GeV^5] // arXiv:1009.5662

  // //  const double sigma02 = 1.881 /constants::hbarc /constants::hbarc; // 18mb - Raju uses 7.2mb
  // const double sigma02 = 0.9/constants::hbarc /constants::hbarc; 
  // const double rt2 =  pow(5.5,2.)/constants::hbarc/constants::hbarc;//(1.69/constants::hbarc/constants::hbarc)*pow(A,2./3.);// assuming R_A = R_0 A^(1/3), with R_0=1.3 fm (R_0^2 = 1.69 fm^2)
  // const double bdep_p = sigma02/2./PI/Bp;
  // const double oomph = 2.21; //for Pb
  // const double bdep_A = oomph*bdep_p;// this is for Pb 
  // //  const double bindep_A = 1./(oomph*(1.69/constants::hbarc/constants::hbarc)/pow(A,1./3.)/2./Bp)*bdep_A; // assuming R_A = R_0 A^(1/3), with R_0=1.3 fm (R_0^2 = 1.69 fm^2)
  // const double bindep_A = A*sigma02/PI/rt2; 
  // const double bdep_fluc_p = bdep_p;
  // const double bdep_fluc_A = bdep_p; // same as for the proton
  const double BKfraction = 0.; //fraction of BK (1- fraction of MV) in BK mode 2 for testing. 
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
  int useFluc;
  int bdep;
  int Ybin;

  //----
  
  double roots;
  double mIR;
  double alphas;
  double Bp;
  double Bq;
  double sigma02;
  double rt2;
  double bdep_p;
  double bdep_A;
  double bindep_A;
  double bdep_fluc_p;
  double bdep_fluc_A;
  double A;
};

// Kinematic variables for ICEM
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
  
// Transverse profiles with fluctuations

double returnTA2D(double x, double y, Glauber *glauberClass, int Ybin){
  return glauberClass->returnNucleusTA(x, y, Ybin);
}

double returnTp2D(double x, double y, Glauber *glauberClass, int Ybin){
  return glauberClass->returnProtonTp(x, y, Ybin);
}

double returnTA(double R, TAInt *TAclass){
return TAclass->returnTA(R);
}

// // Unintegrated gluon distribution for the proton in GBW
// double PhipGBW(double k, double R, double Qs){
//   return constants::CF*k*k*constants::Nc*constants::PI/constants::alphas/constants::CA/Qs/Qs*exp(R*R/(2.*constants::Bp))
//     *exp(-constants::CF*exp(R*R/(2.*constants::Bp))*k*k/(constants::CA*Qs*Qs));
// }

// choose between MV and BK 
double PhipFluc(double k, double Tp, double Qs, double sizeFactor, MV *mv, int BK, double x, double bdep_fluc_p, double alphas){
  if(BK==1){
    Tp = Tp * bdep_fluc_p; // To convert BK b-indepedent to BK b-dependent
    return mv->PhipBKFluc(k, Tp, x, alphas);
  }
  else{
    Tp = Tp * bdep_fluc_p; // To convert MV b-indepedent to MV b-dependent
    return mv->PhipFluc(k, Tp, Qs, sizeFactor, alphas);
  }
}

double Phip(double k, double R, double Qs, double sizeFactor, MV *mv, int BK, double x, int bdep, double bdep_p, double Bp, double alphas){
  double bfactor = 1.; 
  if(bdep==1) // b-dependent
    bfactor = bdep_p; 
 
  if(BK==0){  // MV
    return mv->Phip(k, R, Qs, sizeFactor, bfactor, Bp, alphas);
  }
  else if(BK==1){ // BK
    return mv->PhipBK(k, R, sizeFactor,x, bfactor, Bp, alphas);
  }
  else if (BK==2){
    double rv = constants::BKfraction*mv->PhipBK(k, R, sizeFactor,x, bfactor, Bp, alphas) + (1.-constants::BKfraction)*mv->Phip(k, R, Qs, sizeFactor,bfactor, Bp, alphas);
    return rv;
  }
  else 
    return 0;
}

// // Unintegrated gluon distribution for lead (Not revisited)
// double PhitGBW(double k, double TA, double Qs){
//   return constants::PI*k*k*constants::Nc/constants::alphas*constants::CF*exp(-constants::CF*k*k/(constants::CA*Qs*Qs*TA))/(constants::CA*Qs*Qs*TA);
// }

// choose between BK and MV
double Phit(double k, double TA, double Qs, MV *mv, int BK, double x, int bdep, int useFluc, double bindep_A, double bdep_A, double bdep_fluc_A, double alphas){
  if(BK==1){
      if(useFluc==1){
        TA = TA*bdep_fluc_A;
        return mv->PhitBK(k, TA, x, alphas);
      }
      else{
        if(bdep==1){
            TA = TA*bdep_A;
            return mv->PhitBK(k, TA, x, alphas);
        }
        else{
            TA = TA*bindep_A;  
            return mv->PhitBK(k, TA, x, alphas);
        }
     }
  } 
  else if (BK==0){
    if(useFluc==1){
      TA = TA*bdep_fluc_A;
      return mv->Phit(k, TA, Qs, alphas);
    }
    else{
      if(bdep==1){
        TA = TA*bdep_A;
        return mv->Phit(k, TA, Qs, alphas);
      }
      else{
        TA = TA*bindep_A;  
        return mv->Phit(k, TA, Qs, alphas);
      }
    }
  }
  else if (BK==2){
    if(useFluc==1){
      TA = TA*bdep_fluc_A;
      double rv = constants::BKfraction*mv->PhitBK(k, TA, x, alphas) + (1.-constants::BKfraction)*mv->Phit(k, TA, Qs, alphas);
      return rv;
    }
    else{
      if(bdep==1){
        TA = TA*bdep_A;
        double rv = constants::BKfraction*mv->PhitBK(k, TA, x, alphas) + (1.-constants::BKfraction)*mv->Phit(k, TA, Qs, alphas);
        return rv;
      }
      else{
        TA = TA*bindep_A;  
        double rv = constants::BKfraction*mv->PhitBK(k, TA, x, alphas) + (1.-constants::BKfraction)*mv->Phit(k, TA, Qs, alphas);
        return rv;
      }
    }
  }
  else
    return 0.;
  //return PhitGBW(k, TA, Qs);
}

// FT of fundamental S of the target
double StFGBW(double k, double TA, double Qs){
  return 4.*constants::PI*exp(-(k*k/(Qs*Qs*TA)))/(Qs*Qs*TA);
} 

// choose between MV and BK 
double StF(double k, double TA, double Qs, MV *mv, int BK, double x, int bdep, int useFluc, double bindep_A, double bdep_A, double bdep_fluc_A){
  if(BK==1){
    if(useFluc==1){
      TA = TA*bdep_fluc_A;
      return mv->StFBK(k, TA, x);
    }
    else{
      if(bdep==1){
        TA = TA*bdep_A;
        return mv->StFBK(k, TA, x);
      }
      else{
        TA = TA*bindep_A;
        return mv->StFBK(k, TA, x);
      }
    }
  } 
  else if (BK==0){
    if(useFluc==1){
      TA = TA*bdep_fluc_A;
      return mv->StF(k, TA, Qs);
    }
    else{
      if(bdep==1){
        TA = TA*bdep_A;
        return mv->StF(k, TA, Qs);
      }
      else{
        TA = TA*bindep_A;  
        return mv->StF(k, TA, Qs);
      }
    }
  }
  else if(BK==2){
    if(useFluc==1){
      TA = TA*bdep_fluc_A;
      double rv = constants::BKfraction*mv->StFBK(k, TA, x) + (1.-constants::BKfraction)*mv->StF(k, TA, Qs);
      return rv;
    }
    else{
      if(bdep==1){
        TA = TA*bdep_A;
        double rv = constants::BKfraction*mv->StFBK(k, TA, x) + (1.-constants::BKfraction)*mv->StF(k, TA, Qs);
        return rv;
      }
      else{
        TA = TA*bindep_A;
        double rv = constants::BKfraction*mv->StFBK(k, TA, x) + (1.-constants::BKfraction)*mv->StF(k, TA, Qs);
        return rv;
      }
    }
  }
  else
    return 0.;
}

// NRQCD in what follows

////////////////////////////////////////////////////
///// b-independent J/Psi cross section ///////////
//////////////////////////////////////////////////

static int JPsiIntegrandNRQCDCsNoB(const int *ndim, const cubareal xx[],
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

  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double Bp = static_cast<params*>(userdata)->Bp;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;

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
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);

  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
 
    double myTA = 1.;

    if(pminuskminusk1minuskprime>30.){
        f=0.;
    }
    else{
     f = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
       *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_cs
       *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *sigma02
        *constants::PI*rt2
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

static int JPsiIntegrandNRQCDCoNoB(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobco4k xx[0]
#define nobco4phik xx[1]
#define nobco4k1 xx[2]
#define nobco4phik1 xx[3]
#define nobco4p xx[4]
#define nobco4phip xx[5]

  double kscale = 30.;
  double pscale = 30.;

  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

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
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/roots;
  
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
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
    
    if(pminuskminusk1>30.){
    f=0.;
  }
  else
    {
    double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
    double myTA = 1.;
 
    f = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *sigma02
        *constants::PI*rt2
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
  }
  return 0;
}

//////////////////////////////////////////////////////////
///// b-independent J/Psi cross section for pt spectrum
///// Sum of color singlet + color octet //
/////////////////////////////////////////////////////////

static int JPsiIntegrandNRQCDNoBNoPt(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobpt4k xx[0]
#define nobpt4phik xx[1]
#define nobpt4k1 xx[2]
#define nobpt4phik1 xx[3]
#define nobpt4kprime xx[4]
#define nobpt4phikprime xx[5]
#define nobpt4phip xx[6]


  double kscale = 30.;

  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double p = PT; // don't scale pT 
  double phip = nobpt4phip*2.*constants::PI;
  double k = nobpt4k*kscale;
  double phik = nobpt4phik*2.*constants::PI;
  double k1 = nobpt4k1*kscale;
  double phik1 = nobpt4phik1*2.*constants::PI;
  double kprime = nobpt4kprime*kscale;
  double phikprime = nobpt4phikprime*2.*constants::PI;
  
  double xp = sqrt(4*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
 
  double myTA = 1.;

  double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *sigma02
      *constants::PI*rt2
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 
  
  double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
      +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
      +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
  
  double pminuskminusk1x = px-kx-k1x;
  double pminuskminusk1y = py-ky-k1y;
  double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  
  if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
  }
  else
    {
  double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *sigma02
        *constants::PI*rt2
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 

    f = (singlet + octet)*p*2.*constants::PI;
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
    // d2R
    // d2b
    // pdphi
    // d2k
    // d2k1
    // d2kprime
    }
  }
  return 0;
} 

/////////////////////////////////////////////
///// b-independent J/Psi Mean pt calculation
///// Sum of color singlet + color octet //
//////////////////////////////////////////

// Numerator
static int JPsiIntegrandNRQCDNoBAvPtNum(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobavpt4k xx[0]
#define nobavpt4phik xx[1]
#define nobavpt4k1 xx[2]
#define nobavpt4phik1 xx[3]
#define nobavpt4kprime xx[4]
#define nobavpt4phikprime xx[5]
#define nobavpt4p xx[6]
#define nobavpt4phip xx[7]


  double kscale = 30.;
  double pscale = 30.;

  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double p = nobavpt4p*pscale; 
  double phip = nobavpt4phip*2.*constants::PI;
  double k = nobavpt4k*kscale;
  double phik = nobavpt4phik*2.*constants::PI;
  double k1 = nobavpt4k1*kscale;
  double phik1 = nobavpt4phik1*2.*constants::PI;
  double kprime = nobavpt4kprime*kscale;
  double phikprime = nobavpt4phikprime*2.*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
 
  double myTA = 1.;

  
  double pminuskminusk1x = px-kx-k1x;
  double pminuskminusk1y = py-ky-k1y;
  double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  
  if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
  }
  else
    {
  double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep,bdep_p,Bp, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *sigma02
      *constants::PI*rt2
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 
  
  double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
      +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
      +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);

  double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *sigma02
        *constants::PI*rt2
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 

  f = (singlet + octet)*p*p*pscale*2.*constants::PI;
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
    // d2R
    // d2b
    // pd^2p
    // d2k
    // d2k1
    // d2kprime
    }
  }
  return 0;
} 

// Denominator
static int JPsiIntegrandNRQCDNoBAvPtDen(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobavptden4k xx[0]
#define nobavptden4phik xx[1]
#define nobavptden4k1 xx[2]
#define nobavptden4phik1 xx[3]
#define nobavptden4kprime xx[4]
#define nobavptden4phikprime xx[5]
#define nobavptden4p xx[6]
#define nobavptden4phip xx[7]


  double kscale = 30.;
  double pscale = 30.;

  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;

  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double p = nobavptden4p*pscale; 
  double phip = nobavptden4phip*2.*constants::PI;
  double k = nobavptden4k*kscale;
  double phik = nobavptden4phik*2.*constants::PI;
  double k1 = nobavptden4k1*kscale;
  double phik1 = nobavptden4phik1*2.*constants::PI;
  double kprime = nobavptden4kprime*kscale;
  double phikprime = nobavptden4phikprime*2.*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
 
  double myTA = 1.;
  
  double pminuskminusk1x = px-kx-k1x;
  double pminuskminusk1y = py-ky-k1y;
  double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  
  double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
      +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
      +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
  
  
  if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
  }
  else
    {

  double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *sigma02
      *constants::PI*rt2
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 

  double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, 0., Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *sigma02
        *constants::PI*rt2
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 

  f = (singlet + octet)*p*pscale*2.*constants::PI;
    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
    // d2R
    // d2b
    // pd^2p
    // d2k
    // d2k1
    // d2kprime
    }
  }
  return 0;
} 
//////////////////////////////////////////////////
///// b-dependent J/Psi cross section ///////////
/////////////////////////////////////////////////

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

  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

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
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
 
  double myTA = returnTA(Rminusb,TAclass); 

  if(pminuskminusk1minuskprime>30.){
    f=0.;
  }
  else{
    f = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
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

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 2./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!
  
  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

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
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else{
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
    
    double myTA = returnTA(Rminusb,TAclass); 
    if(pminuskminusk1>30.){
    f=0.;
  }
  else
    {
    
    f = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
      *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
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
  }
  return 0;
} 

////////////////////////////////////////////////////////
///// b-dependent J/Psi cross section for pt spectrum //
///// Sum of color singlet + color octet //
///////////////////////////////////////////////////////

static int JPsiIntegrandNRQCDNoPT(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define noptR xx[0]
#define noptphiR xx[1]
#define noptb xx[2]
#define noptphib xx[3]
#define nopt4k xx[4]
#define nopt4phik xx[5]
#define nopt4k1 xx[6]
#define nopt4phik1 xx[7]
#define nopt4kprime xx[8]
#define nopt4phikprime xx[9]
#define nopt4phip xx[10]

  double kscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!

  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double R = noptR*Rscale;
  double b = noptb*bscale;
  double p = PT;
  double phip = nopt4phip*2.*constants::PI;
  double k = nopt4k*kscale;
  double phik = nopt4phik*2.*constants::PI;
  double k1 = nopt4k1*kscale;
  double phik1 = nopt4phik1*2.*constants::PI;
  double kprime = nopt4kprime*kscale;
  double phikprime = nopt4phikprime*2.*constants::PI;
  double phiR = noptphiR*2*constants::PI;
  double phib = noptphib*2*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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

    double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
  
    double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);
 
    double myTA = returnTA(Rminusb,TAclass); //(2.37=0.4*(208)^(1/3))
  
  
    double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_cs
        *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *R*Rscale*2.*constants::PI
        *b*bscale*2.*constants::PI
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI
        *kprime*kscale*2.*constants::PI; 

    double pminuskminusk1x = px-kx-k1x;
    double pminuskminusk1y = py-ky-k1y;
    double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  
    double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
    if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
  }
  else
    {
    double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
       *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *R*Rscale*2.*constants::PI
        *b*bscale*2.*constants::PI
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI;

    f = (singlet + octet)*p*2.*constants::PI;

    // scaled momenta above (in PT)
    // last rows are scaling of integration measures:
    // d2R
    // d2b
    // pdphi
    // d2k
    // d2k1
    // d2kprime
    }
   }
  return 0;
} 

/////////////////////////////////////////////
///// b-dependent J/Psi Mean pt calculation/
///// Sum of color singlet + color octet //
//////////////////////////////////////////

// Numerator
static int JPsiIntegrandNRQCDAvPtNum(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define avptnumR xx[0]
#define avptnumphiR xx[1]
#define avptnumb xx[2]
#define avptnumphib xx[3]
#define avptnum4k xx[4]
#define avptnum4phik xx[5]
#define avptnum4k1 xx[6]
#define avptnum4phik1 xx[7]
#define avptnum4phip xx[8]
#define avptnum4p xx[9]
#define avptnum4kprime xx[10]
#define avptnum4phikprime xx[11]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 2./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!
  
  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double R = avptnumR*Rscale;
  double b = avptnumb*bscale;
  double p = avptnum4p*pscale;
  double phip = avptnum4phip*2.*constants::PI;
  double k = avptnum4k*kscale;
  double phik = avptnum4phik*2.*constants::PI;
  double k1 = avptnum4k1*kscale;
  double phik1 = avptnum4phik*2.*constants::PI;
  double kprime = avptnum4kprime*kscale;
  double phikprime = avptnum4phikprime*2.*constants::PI;
  double phiR = avptnumphiR*2*constants::PI;
  double phib = avptnumphib*2*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f = 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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

    double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
  
    double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
    double myTA = returnTA(Rminusb,TAclass); 


    double kprimex = kprime*cos(phikprime);
    double kprimey = kprime*sin(phikprime);
    double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
    double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
    double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                    +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
    double phi_pminuskminusk1minuskprime = atan2(pminuskminusk1minuskprimey,pminuskminusk1minuskprimex);

    if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
    }
    else
    {

    double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);

    double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep,bdep_p,Bp, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *sigma02
      *constants::PI*rt2
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 

  
    double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *R*Rscale*2.*constants::PI
        *b*bscale*2.*constants::PI
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 


    f = (singlet + octet)*p*p*pscale*2.*constants::PI;
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // pd2p
  // d2k
  // d2k1
    }
   
  }
  return 0;
} 

// Denominator
static int JPsiIntegrandNRQCDAvPtDen(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define avptdenR xx[0]
#define avptdenphiR xx[1]
#define avptdenb xx[2]
#define avptdenphib xx[3]
#define avptden4k xx[4]
#define avptden4phik xx[5]
#define avptden4k1 xx[6]
#define avptden4phik1 xx[7]
#define avptden4phip xx[8]
#define avptden4p xx[9]
#define avptden4kprime xx[10]
#define avptden4phikprime xx[11]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 2./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; //!!
  
  int BK = static_cast<params*>(userdata)->BK;
  double Y = static_cast<params*>(userdata)->Y;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double R = avptdenR*Rscale;
  double b = avptdenb*bscale;
  double p = avptden4p*pscale;
  double phip = avptden4phip*2.*constants::PI;
  double k = avptden4k*kscale;
  double phik = avptden4phik*2.*constants::PI;
  double k1 = avptden4k1*kscale;
  double phik1 = avptden4phik*2.*constants::PI;
  double kprime = avptden4kprime*kscale;
  double phikprime = avptden4phikprime*2.*constants::PI;
  double phiR = avptdenphiR*2*constants::PI;
  double phib = avptdenphib*2*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f = 0;
  }
  else{
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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

    double Rminusb = sqrt(R*R+b*b-2.*R*b*cos(phiR-phib));
  
    double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
    double myTA = returnTA(Rminusb,TAclass); 


    double kprimex = kprime*cos(phikprime);
    double kprimey = kprime*sin(phikprime);
    double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
    double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
    double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                    +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
    double phi_pminuskminusk1minuskprime = atan2(pminuskminusk1minuskprimey,pminuskminusk1minuskprimex);

    double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);
    
  if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
  }
  else
    {

    double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep,bdep_p,Bp, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(kprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(pminuskminusk1minuskprime,myTA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *sigma02
      *constants::PI*rt2
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 

  
    double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *Phip(k1, R, Qsp, sizeFactor,mv,BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *R*Rscale*2.*constants::PI
        *b*bscale*2.*constants::PI
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 


    f = (singlet + octet)*p*pscale*2.*constants::PI;
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // pd2p
  // d2k
  // d2k1
    }
   
  }
  return 0;
} 

//////////////////////////////////////////////////
///// Fluctuating b J/Psi cross section //////////
/////////////////////////////////////////////////

//////////////////////////////////////////////////
///// Fluctuating b J/Psi cross section //////////
/////////////////////////////////////////////////

static int JPsiIntegrandNRQCDAvPtNumFluc(const int *ndim, const cubareal xx[],
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

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;
  
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
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    // get sums of vectors
    double px = p*cos(phip); 
    double py = p*sin(phip);
    double kx = k*cos(phik);
    double ky = k*sin(phik);
    double k1x = k1*cos(phik1);
    double k1y = k1*sin(phik1);
    double kprimex = kprime*cos(phikprime);
    double kprimey = kprime*sin(phikprime);

    double pminuskminusk1x = px-kx-k1x;
    double pminuskminusk1y = py-ky-k1y;
    double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);

    
    double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
    double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
    double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                            +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
    
    double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);

    double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
    
    double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
    double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
    

  if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
  }
  else
    {

    double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *PhipFluc(k1, Tp, Qsp, sizeFactor,mv,BK,xp,bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(kprime,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(pminuskminusk1minuskprime,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *Rscale*Rscale
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 
      
  
    double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *PhipFluc(k1, Tp, Qsp, sizeFactor,mv,BK,xp,bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *Rscale*Rscale
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 

    
 
    f = (singlet + octet)*p*p*pscale*2.*constants::PI;

    }

  }
  return 0;
}

static int JPsiIntegrandNRQCDAvPtDenFluc(const int *ndim, const cubareal xx[],
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

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

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
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    // get sums of vectors
    double px = p*cos(phip); 
    double py = p*sin(phip);
    double kx = k*cos(phik);
    double ky = k*sin(phik);
    double k1x = k1*cos(phik1);
    double k1y = k1*sin(phik1);
    double kprimex = kprime*cos(phikprime);
    double kprimey = kprime*sin(phikprime);

    double pminuskminusk1x = px-kx-k1x;
    double pminuskminusk1y = py-ky-k1y;
    double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);

    double pminuskminusk1minuskprimex = px-kx-k1x-kprimex;
    double pminuskminusk1minuskprimey = py-ky-k1y-kprimey;
    double pminuskminusk1minuskprime = sqrt(pminuskminusk1minuskprimex*pminuskminusk1minuskprimex
                                            +pminuskminusk1minuskprimey*pminuskminusk1minuskprimey);
    
    double H_cs = constants::ldme_singlet*nrqcd::singlet(p, phip, k1, phik1,kprime, phikprime, k, phik,m);
    
    double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);

    double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
    double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
    
  if(pminuskminusk1>30. || pminuskminusk1minuskprime > 30.){
    f=0.;
  }
  else
    {

    double singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *PhipFluc(k1, Tp, Qsp, sizeFactor,mv,BK,xp,bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_cs
      *(StF(k,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(kprime,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA
      *StF(pminuskminusk1minuskprime,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *Rscale*Rscale
      *k*kscale*2.*constants::PI
      *k1*kscale*2.*constants::PI
      *kprime*kscale*2.*constants::PI; 
      
  
    double octet = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *PhipFluc(k1, Tp, Qsp, sizeFactor,mv,BK,xp,bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *Rscale*Rscale
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 

    f = (singlet + octet)*p*pscale*2.*constants::PI;

    }

  }
  return 0;
}

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

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

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
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
    
    double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
    double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
    
    if(pminuskminusk1minuskprime>30.){
      f=0.;
    }
    else
      {
        f = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
          *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_cs
          *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(kprime,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1minuskprime,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
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

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double roots = static_cast<params*>(userdata)->roots;
  double alphas = static_cast<params*>(userdata)->alphas;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

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
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
 
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
 
  if(pminuskminusk1>30.){
    f=0.;
  }
  else
    {
      f = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
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


///////////////////////////////////////////////////////////
///// Fluctuating b J/Psi cross section for pt spectrum ///
//////////////////////////////////////////////////////////

static int JPsiIntegrandNRQCDFlucNoPT(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define noptfRx xx[0]
#define noptfRy xx[1]
#define noptf4k xx[2]
#define noptf4phik xx[3]
#define noptf4k1 xx[4]
#define noptf4phik1 xx[5]
#define noptf4kprime xx[6]
#define noptf4phikprime xx[7]
#define noptf4phip xx[8]

  double kscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;  
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double roots = static_cast<params*>(userdata)->roots;
  double alphas = static_cast<params*>(userdata)->alphas;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double Rx = noptfRx*Rscale-Rscale/2.;
  double Ry = noptfRy*Rscale-Rscale/2.;
  double p = PT;
  double phip = noptf4phip*2.*constants::PI;
  double k = noptf4k*kscale;
  double phik = noptf4phik*2.*constants::PI;
  double k1 = noptf4k1*kscale;
  double phik1 = noptf4phik1*2.*constants::PI;
  double kprime = noptf4kprime*kscale;
  double phikprime = noptf4phikprime*2.*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
    
    double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
    double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);

    
    double singlet =0.;
    if(pminuskminusk1minuskprime>30.){
      singlet=0.;
    }
    else
      {
        singlet = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
          *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_cs
          *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(kprime,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1minuskprime,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
          *Rscale*Rscale
          *k*kscale*2.*constants::PI
          *k1*kscale*2.*constants::PI
          *kprime*kscale*2.*constants::PI;
      }
    
    double pminuskminusk1x = px-kx-k1x;
    double pminuskminusk1y = py-ky-k1y;
    double pminuskminusk1 = sqrt(pminuskminusk1x*pminuskminusk1x
                                    +pminuskminusk1y*pminuskminusk1y);
  
    double H_co = constants::ldme_octet_s10*nrqcd::octets10(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_s13*nrqcd::octets13(p, phip, k1, phik1, k, phik,m)
                +constants::ldme_octet_p3j*nrqcd::octetp3j(p, phip, k1, phik1, k, phik,m);
 
    double octet=0.;
    if(pminuskminusk1>30.){
      octet=0.;
    }
    else
      {
        octet =  alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
          *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_co
          *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
          *Rscale*Rscale
          *k*kscale*2.*constants::PI
          *k1*kscale*2.*constants::PI;
      }
    
    f = (singlet + octet)*p*2.*constants::PI;
  
  }
  return 0;
} 

static int JPsiIntegrandNRQCDCsFlucNoPT(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define noptfcsRx xx[0]
#define noptfcsRy xx[1]
#define noptfcs4k xx[2]
#define noptfcs4phik xx[3]
#define noptfcs4k1 xx[4]
#define noptfcs4phik1 xx[5]
#define noptfcs4kprime xx[6]
#define noptfcs4phikprime xx[7]
#define noptfcs4phip xx[8]

  double kscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;  
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double roots = static_cast<params*>(userdata)->roots;
  double alphas = static_cast<params*>(userdata)->alphas;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double Rx = noptfcsRx*Rscale-Rscale/2.;
  double Ry = noptfcsRy*Rscale-Rscale/2.;
  double p = PT;
  double phip = noptfcs4phip*2.*constants::PI;
  double k = noptfcs4k*kscale;
  double phik = noptfcs4phik*2.*constants::PI;
  double k1 = noptfcs4k1*kscale;
  double phik1 = noptfcs4phik1*2.*constants::PI;
  double kprime = noptfcs4kprime*kscale;
  double phikprime = noptfcs4phikprime*2.*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
    
    double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
    double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
    
    if(pminuskminusk1minuskprime>30.){
      f=0.;
    }
    else
      {
        f = alphas/(pow(2.*constants::PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
          *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_cs
          *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(kprime,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1minuskprime,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
          *Rscale*Rscale
          *p*2.*constants::PI
          *k*kscale*2.*constants::PI
          *k1*kscale*2.*constants::PI
          *kprime*kscale*2.*constants::PI; 
        // scaled momenta above (in PT)
        // last rows are scaling of integration measures:
        // dRxdRy
        // 2pi p
        // d2k
        // d2k1
        // d2kprime
      }
  }
  return 0;
} 


static int JPsiIntegrandNRQCDCoFlucNoPT(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define noptfcoRx xx[0]
#define noptfcoRy xx[1]
#define noptfco4k xx[2]
#define noptfco4phik xx[3]
#define noptfco4k1 xx[4]
#define noptfco4phik1 xx[5]
#define noptfco4phip xx[6]

  double kscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables  
  double Rx = noptfcoRx*Rscale-Rscale/2.;
  double Ry = noptfcoRy*Rscale-Rscale/2.;
  double p = PT;
  double phip = noptfco4phip*2.*constants::PI;
  double k = noptfco4k*kscale;
  double phik = noptfco4phik*2.*constants::PI;
  double k1 = noptfco4k1*kscale;
  double phik1 = noptfco4phik1*2.*constants::PI;
  
  double xp = sqrt(4.*m*m+p*p)*exp(Y)/roots;
  double xA = sqrt(4.*m*m+p*p)*exp(-Y)/roots;
  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else {
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
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
 
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
 
  if(pminuskminusk1>30.){
    f=0.;
  }
  else
    {
      f = alphas/(pow(2.*constants::PI,7.)*(double(constants::Nc)*double(constants::Nc)-1.))
        *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H_co
        *(StF(k,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(pminuskminusk1,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
        *Rscale*Rscale
        *p*2.*constants::PI
        *k*kscale*2.*constants::PI
        *k1*kscale*2.*constants::PI; 

  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // dRxdRy
  // 2pi p
  // d2k
  // d2k1
    }
  }
  return 0;
}

///////////////////////////////////////////////////
///// b-independent J/Psi cross section ///////////
///////////////////////////////////////////////////

static int JPsiIntegrandAllNoB(const int *ndim, const cubareal xx[],
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
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

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
  in.PT  = PT*M/constants::mJPsi;
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

  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/roots;

  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else{
    double betax = (M/constants::mJPsi*PT)/sqrt(M*M+(M/constants::mJPsi*PT)*(M/constants::mJPsi*PT));
    //double betax = PT/sqrt(M*M+PT*PT);
    double gammax = 1./sqrt(1.-betax*betax);
    // get Jacobian  
    double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));
    
    double myTA = 1.;
    double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);

    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

    f = alphas*double(constants::Nc)*double(constants::Nc)
      /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, 0., Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp/(k1*k1)*H*J
      *(StF(pplusqminusk1minusk,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(k,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
      *sigma02
      *constants::PI*rt2
      //*2.*constants::PI*constants::Bp
      //*2.*constants::PI*constants::Bt
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

//////////////////////////////////////////////////
///// b-dependent J/Psi cross section ///////////
/////////////////////////////////////////////////

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


  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 12./constants::hbarc; 
  // Qs will be made rapidity dependent
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + qqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = qqqtilde*qtildescale;
  double PT = qqPT*pscale; 
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
  in.PT  = PT*M/constants::mJPsi;
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
  
  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/roots;

  
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }

  else if (xA>1.){
    f= 0;
  }
  else{
    //    double betax = PT/sqrt(M*M+PT*PT);
    double betax = (M/constants::mJPsi*PT)/sqrt(M*M+(M/constants::mJPsi*PT)*(M/constants::mJPsi*PT));

    double gammax = 1./sqrt(1.-betax*betax);
    // get Jacobian  
    double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));
    double myTA = returnTA(Rminusb,TAclass); 
    double H = Hard::all(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, k, phik, yp, yq, m);

    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

    f = alphas*double(constants::Nc)*double(constants::Nc)
      /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
      *Phip(k1, R, Qsp, sizeFactor, mv, BK,xp,bdep, bdep_p, Bp, alphas)*factorxp/(k1*k1)*H*J
      *(StF(pplusqminusk1minusk,myTA,QsA,mv, BK, xA, bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(k,myTA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA)
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

//////////////////////////////////////////////////
///// Fluctuating b J/Psi cross section //////////
/////////////////////////////////////////////////

static int JPsiIntegrandICEMNumFluc(const int *ndim, const cubareal xx[],
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
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
 
  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
 
  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + fqqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = fqqqtilde*qtildescale;
  double PT = fqqPT*pscale; 

  double Rx = fqqRx*Rscale-Rscale/2.;
  double Ry = fqqRy*Rscale-Rscale/2.;
  double k = fqq4k*kscale;
  double phik = fqq4phik*2.*constants::PI;
  double k1 = fqq4k1*kscale;
  double phik1 = fqq4phik1*2.*constants::PI;

  kinPair in;
  in.M = M;
  //  in.PT  = PT;
  in.PT  = PT*(M/constants::mJPsi); // evaluate cross section at PT' =  PT*(M/constants::mJPsi) (see (6) in https://arxiv.org/pdf/1609.06042.pdf)
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

  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/roots;


  if (xp>1.){
        f = 0.;
  }
  else if (xA>1.){
       f = 0.;
  }
  
  else{
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);

  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

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

  // get Jacobian
  //  double betax = PT/sqrt(M*M+PT*PT);
  double betax = (M/constants::mJPsi*PT)/sqrt(M*M+(M/constants::mJPsi*PT)*(M/constants::mJPsi*PT));
  double gammax = 1./sqrt(1.-betax*betax);
  
  double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));

  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);

  // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
  f = alphas*double(constants::Nc)*double(constants::Nc)
    /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H*J
    *(StF(pplusqminusk1minusk,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(k,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A))*factorxA
    *Rscale*Rscale
    *PT*PT*pscale*2.*constants::PI
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
  } 
  return 0;
}  

static int JPsiIntegrandICEMDenFluc(const int *ndim, const cubareal xx[],
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
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
 
  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
 
  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + fqqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = fqqqtilde*qtildescale;
  double PT = fqqPT*pscale; 

  double Rx = fqqRx*Rscale-Rscale/2.;
  double Ry = fqqRy*Rscale-Rscale/2.;
  double k = fqq4k*kscale;
  double phik = fqq4phik*2.*constants::PI;
  double k1 = fqq4k1*kscale;
  double phik1 = fqq4phik1*2.*constants::PI;

  kinPair in;
  in.M = M;
  in.PT  = PT*(M/constants::mJPsi); // evaluate cross section at PT' =  PT*(M/constants::mJPsi) (see (6) in https://arxiv.org/pdf/1609.06042.pdf)
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

  double xp = (sqrt(p*p+m*m)*exp(yp)+sqrt(q*q+m*m)*exp(yq))/roots;
  double xA = (sqrt(p*p+m*m)*exp(-yp)+sqrt(q*q+m*m)*exp(-yq))/roots;


  if (xp>1.){
        f = 0.;
  }
  else if (xA>1.){
       f = 0.;
  }
  
  else{
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);

  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

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

  // get Jacobian
  //double betax = PT/sqrt(M*M+PT*PT);
  double betax = (M/constants::mJPsi*PT)/sqrt(M*M+(M/constants::mJPsi*PT)*(M/constants::mJPsi*PT));
  double gammax = 1./sqrt(1.-betax*betax);
  
  double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));

  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);

  // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
  f = alphas*double(constants::Nc)*double(constants::Nc)
    /(2.*pow(2.*constants::PI,10.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *PhipFluc(k1, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp/(k1*k1)*H*J
    *(StF(pplusqminusk1minusk,TA,QsA,mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A)*factorxA*StF(k,TA,QsA,mv,BK,xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A))*factorxA
    *Rscale*Rscale
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
  } 
  return 0;
}  


// Gluons

///////////////////////////////////////////////////
///// b-independent gluons cross section //////////
//////////////////////////////////////////////////

static int GluonsNoB(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobk xx[0]
#define nobphik xx[1]
#define nobp xx[2]
#define nobphi xx[3]

  double kscale = 30.;
  double pscale = 30.;

  double lambda = static_cast<params*>(userdata)->lambda;
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
    
  double k = nobk*kscale;
  double p = nobp*pscale+lambda;

  double phi = 2.*constants::PI*nobphi;
  double phik = 2.*constants::PI*nobphik;

  double xp = (nobp*pscale+lambda)*exp(Y)/roots;
  double xA = (nobp*pscale+lambda)*exp(-Y)/roots;
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

  double TA = 1.; // To avoid impact parameter dependence. We also set R=0 inside Phip for the same purpose
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    f = alphas/constants::CF/p/p/pow((2*constants::PI*constants::PI),3.)
      *Phip(k, 0, Qsp, sizeFactor, mv, BK, xp,bdep, bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt(p*p + k*k - 2.*p*k*cos(phi-phik)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)*factorxA
      *2.*constants::PI*k*kscale  //kdkdphik
      *sigma02  //R-integral
      *constants::PI*rt2  // b-integral
      *2.*constants::PI*pscale*p; //pdpdphip
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
    return 1;
}


static int GluonsNoBNoPt(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define nobnoptk xx[0]
#define nobnoptphik xx[1]
#define nobnoptphi xx[2]

  double kscale = 30.;
  double pscale = 30.;

  double lambda = static_cast<params*>(userdata)->lambda;
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double PT = static_cast<params*>(userdata)->PT;
    
  double k = nobnoptk*kscale;
  double p = PT;

  double phi = 2.*constants::PI*nobnoptphi;
  double phik = 2.*constants::PI*nobnoptphik;

  double xp = (p)*exp(Y)/roots;
  double xA = (p)*exp(-Y)/roots;
  double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
  double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

  double TA = 1.; // To avoid impact parameter dependence. We also set R=0 inside Phip for the same purpose
  double factorxp = pow(1.-xp,4.);
  double factorxA = pow(1.-xA,4.);
  if (xp>1.){
    f = 0;
  }
  else if (xA>1.){
    f= 0;
  }
  else {
    f = alphas/constants::CF/(p*p+mIR*mIR)/pow((2*constants::PI*constants::PI),3.)
      *Phip(k, 0, Qsp, sizeFactor, mv, BK, xp,bdep, bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt(p*p + k*k - 2.*p*k*cos(phi-phik)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)*factorxA
      *2.*constants::PI*k*kscale  //kdkdphik
      *sigma02  //R-integral
      *constants::PI*rt2  // b-integral
      *2.*constants::PI*p; //p dphip
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
    return 1;
}

//////////////////////////////////////////////////
///// b-dependent gluons cross section ///////////
//////////////////////////////////////////////////

static int Gluons(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define gk xx[0]
#define gphik xx[1]
#define gR xx[2]
#define gphiR xx[3]
#define gb xx[4]
#define gphib xx[5]
#define gphi xx[6]
#define gp xx[7]

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 2./constants::hbarc;
  double bscale = 12./constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
 
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*gphi;
  double phik = 2.*constants::PI*gphik;
  double phiR = 2.*constants::PI*gphiR;
  double phib = 2.*constants::PI*gphib;
  double k = gk*kscale;
  double p = gp*pscale+lambda;
  double R = gR*Rscale;
  double b = gb*bscale;

  double xp = (p)*exp(Y)/roots;
  double xA = (p)*exp(-Y)/roots;
  double factorxA = pow(1.-xA,4.);
  double factorxp = pow(1.-xp,4.);

  if (xp>1.){
    f = 0.;
  }
  else if (xA>1.){
    f = 0.;
  }
  else{
   double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
   double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

   double TA = returnTA(sqrt(max(R*R + b*b - 2.*R*b*cos((phiR - phib)),0.)),TAclass);

   f = alphas/constants::CF/(p*p+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
     *Phip(k, R, Qsp, sizeFactor, mv, BK,xp,bdep, bdep_p,Bp, alphas)
     *Phit(sqrt((p)*(p) + k*k - 2.*p*k*cos((phi - phik))), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
       *2.*constants::PI*kscale*k  //kdkdphik
       *2.*constants::PI*Rscale*R  //RdRdphiR
       *2.*constants::PI*bscale*b  //bdbdphib
       *2.*constants::PI*pscale*p //pdpdphip
       *factorxp*factorxA;
   //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}


static int GluonsNoPt(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define gk xx[0]
#define gphik xx[1]
#define gR xx[2]
#define gphiR xx[3]
#define gb xx[4]
#define gphib xx[5]
#define gphi xx[6]

  double kscale = 30.;
  double Rscale = 2./constants::hbarc;
  double bscale = 12./constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double PT = static_cast<params*>(userdata)->PT;
 
  double xp = (PT)*exp(Y)/roots;
  double xA = (PT)*exp(-Y)/roots;
  double factorxA = pow(1.-xA,4.);
  double factorxp = pow(1.-xp,4.);
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*gphi;
  double phik = 2.*constants::PI*gphik;
  double phiR = 2.*constants::PI*gphiR;
  double phib = 2.*constants::PI*gphib;
  double k = gk*kscale;
  double p = PT;
  double R = gR*Rscale;
  double b = gb*bscale;

  if (xp>1.){
    f = 0.;
  }
  else if (xA>1.){
    f = 0.;
  }
  else{
   double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
   double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);

   double TA = returnTA(sqrt(max(R*R + b*b - 2.*R*b*cos((phiR - phib)),0.)),TAclass);

   // if(k>p) //!! cutoff
   //   f=0.;
   // else
   f = alphas/constants::CF/(p*p+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.)//!!
       *Phip(k, R, Qsp, sizeFactor, mv, BK,xp,bdep, bdep_p,Bp, alphas)
       *Phit(sqrt(p*p + k*k - 2.*p*k*cos((phi - phik))), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
       *2.*constants::PI*kscale*k  //kdkdphik
       *2.*constants::PI*Rscale*R  //RdRdphiR
       *2.*constants::PI*bscale*b  //bdbdphib
       *2.*constants::PI*p //pdpdphip
       *factorxp*factorxA;
   //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}

////////////////////////////////////////////////////
///// Fluctuating b gluons cross section ///////////
///////////////////////////////////////////////////

static int GluonsFluc(const int *ndim, const cubareal xx[],
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
  double Rx = fgRx*Rscale-Rscale/2.;
  double Ry = fgRy*Rscale-Rscale/2.;

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;
  double lambda = static_cast<params*>(userdata)->lambda;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double Y = static_cast<params*>(userdata)->Y;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double R = sqrt(Rx*Rx+Ry*Ry);
  double b = sqrt(bx*bx+by*by);
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);

  double p = fgp*pscale+lambda;
  double xp = p*exp(Y)/roots;
  double xA = p*exp(-Y)/roots;
  
  double factorxA = pow(1.-xA,4.);
  double factorxp = pow(1.-xp,4.);
  if (xp>1.){
    f = 0.;
  }
  else if (xA>1.){
    f = 0.;
  }
  else{
    
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
    f = alphas/constants::CF/(p*p+mIR*mIR)/pow((2*constants::PI*constants::PI),3.)
      *PhipFluc(fgk*kscale, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp
      *Phit(sqrt((fgp*pscale+lambda)*(fgp*pscale+lambda) + fgk*fgk*kscale*kscale - 2.*(fgp*pscale+lambda)*fgk*kscale*cos((fgphi - fgphik)*2.*constants::PI)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)*factorxA
      *2.*constants::PI*fgk*kscale*kscale  //kdkdphik
      *Rscale*Rscale  //dRxdRy
      //    *bscale*bscale  //bdxdby
      *2.*constants::PI*pscale*(fgp*pscale+lambda); //pdpdphip
  }
  
  return 1;
}
/////////////////////////////////////////////////////
//////// b independent hadron cross section ///////////
/////////////////////////////////////////////////////
static int HadronsNoB(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define hnobk xx[0]
#define hnobphik xx[1]
#define hnobp xx[2]
#define hnobphi xx[3]

  double z = xx[4];
  double mh = 0.2; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;  
  double mIR = static_cast<params*>(userdata)->mIR;
  
  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double A = static_cast<params*>(userdata)->A;


  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*hnobphi;
  double phik = 2.*constants::PI*hnobphik;

  double k = hnobk*kscale;
  double p = hnobp*pscale+lambda;

  double eta = Y;
  double pg = p/z;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  //assuming yg=yh

  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;

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
    
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    double TA = 1; // To avoid impact parameter dependence. We also set R=0 inside Phip for the same purpose
    
    f = Dh*J* 
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      *Phip(k, 0, Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos(phi - phik)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
      *factorxA
      *2.*constants::PI*kscale*k  //kdkdphik
      *2.*constants::PI*pscale*p //pdpdphip
      *sigma02  //R-integral
      *constants::PI*rt2;  // b-integral
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}

/////////////////////////////////////////////////////
//////// b independent hadron pt spectrum ///////////
/////////////////////////////////////////////////////

static int HadronsNoBNoPt(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define hnobnoptk xx[0]
#define hnobnoptphik xx[1]
#define hnobnoptphi xx[2]

  double z = xx[3];
  double mh = 0.2; //GeV

  double kscale = 30.;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  double mIR = static_cast<params*>(userdata)->mIR;
  
  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double A = static_cast<params*>(userdata)->A;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*hnobnoptphi;
  double phik = 2.*constants::PI*hnobnoptphik;

  double k = hnobnoptk*kscale;
  double p = PT;

  double eta = Y;
  double pg = p/z;
  

  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  // //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  // double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;

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
    
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    double TA = 1; // To avoid impact parameter dependence. We also set R=0 inside Phip for the same purpose
    
    f = Dh*J* 
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.)
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      *Phip(k, 0, Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos(phi - phik)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
      *factorxA
      *2.*constants::PI*kscale*k  //kdkdphik
      *2.*constants::PI*p //pdphip
      *sigma02  //R-integral
      *constants::PI*rt2;  // b-integral
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}

/////////////////////////////////////////////////////
//////// b independent hadron <pt> ///////////
/////////////////////////////////////////////////////

// Numerator 

static int HadronsNoBAvPtNum(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define hnobavptnumk xx[0]
#define hnobavptnumphik xx[1]
#define hnobavptnump xx[2]
#define hnobavptnumphi xx[3]

  double z = xx[4];
  double mh = 0.2; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double mIR = static_cast<params*>(userdata)->mIR;
  
  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double A = static_cast<params*>(userdata)->A;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*hnobavptnumphi;
  double phik = 2.*constants::PI*hnobavptnumphik;

  double k = hnobavptnumk*kscale;
  double p = hnobavptnump*pscale+lambda;

  double eta = Y;
  double pg = p/z;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  // //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  // double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;

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
    
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    double TA = 1; // To avoid impact parameter dependence. We also set R=0 inside Phip for the same purpose
    
    f = Dh*J*
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.)
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.)
      *Phip(k, 0, Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos(phi - phik)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
      *factorxA
      *2.*constants::PI*kscale*k  //kdkdphik
      *2.*constants::PI*pscale*p*p //p^2dpdphip
      *sigma02  //R-integral
      *constants::PI*rt2;  // b-integral
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}


/////////////////////////////////////////////////////
//////// b dependent hadron cross section ///////////
/////////////////////////////////////////////////////

static int Hadrons(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define hk xx[0]
#define hphik xx[1]
#define hR xx[2]
#define hphiR xx[3]
#define hb xx[4]
#define hphib xx[5]
#define hphi xx[6]
#define hp xx[7]

  double z = xx[8];
  double mh = 0.2; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 2./constants::hbarc;
  double bscale = 12./constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;  
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double A = static_cast<params*>(userdata)->A;
 
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*hphi;
  double phik = 2.*constants::PI*hphik;
  double phiR = 2.*constants::PI*hphiR;
  double phib = 2.*constants::PI*hphib;
  double k = hk*kscale;
  double p = hp*pscale+lambda;
  double R = hR*Rscale;
  double b = hb*bscale;

  double eta = Y;
  double pg = p/z;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  // //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  // double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;

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
    
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    double TA = returnTA(sqrt(max(R*R + b*b - 2.*R*b*cos((phiR - phib)),0.)),TAclass);
    
    f = Dh*J*
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.)//!!
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      *Phip(k, R, Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt((pg)*(pg) + k*k - 2.*pg*k*cos((phi - phik))), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
      *factorxA
      *2.*constants::PI*kscale*k  //kdkdphik
      *2.*constants::PI*Rscale*R  //RdRdphiR
      *2.*constants::PI*bscale*b  //bdbdphib
      *2.*constants::PI*pscale*p; //pdpdphip
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}

/////////////////////////////////////////////////////
//////// b dependent hadron pt spectrum ///////////
/////////////////////////////////////////////////////

static int HadronsNoPt(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define hnoptk xx[0]
#define hnoptphik xx[1]
#define hnoptR xx[2]
#define hnoptphiR xx[3]
#define hnoptb xx[4]
#define hnoptphib xx[5]
#define hnoptphi xx[6]

  double z = xx[7];
  
  double mh = 0.2; //GeV

  double kscale = 30.;
  double Rscale = 2./constants::hbarc;
  double bscale = 12./constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double PT = static_cast<params*>(userdata)->PT;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double A = static_cast<params*>(userdata)->A;
 
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*hnoptphi;
  double phik = 2.*constants::PI*hnoptphik;
  double phiR = 2.*constants::PI*hnoptphiR;
  double phib = 2.*constants::PI*hnoptphib;
  double k = hnoptk*kscale;
  double p = PT;
  double R = hnoptR*Rscale;
  double b = hnoptb*bscale;

  double eta = Y;
  double pg = p/z;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);

  // //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO   
  // //  double Dh = kkp::KKPFragmentation(7, 1, z, pg, gluon);
  // double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);

  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;

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
    
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    double TA = returnTA(sqrt(max(R*R + b*b - 2.*R*b*cos((phiR - phib)),0.)),TAclass);
  
   // if(k>p) //!! cutoff
   //   f=0.;
   // else
    // f = Dh/z/z*J*
    //   alphas/constants::CF/(pg*pg+0.4*0.4)/pow((2.*constants::PI*constants::PI),3.) //!!
    //   *Phip(k, R, Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp
    //   *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos((phi - phik))), TA, QsA, mv, BK, xA, bdep, useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
    //   *factorxA
    //   *2.*constants::PI*kscale*k  //kdkdphik
    //   *2.*constants::PI*Rscale*R  //RdRdphiR
    //   *2.*constants::PI*bscale*b  //bdbdphib
    //   *2.*constants::PI*p; //pdphip

    f = Dh*J*
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      *Phip(k, R, Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos((phi - phik))), TA, QsA, mv, BK, xA, bdep, useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
      *factorxA
      *2.*constants::PI*kscale*k  //kdkdphik
      *2.*constants::PI*Rscale*R  //RdRdphiR
      *2.*constants::PI*bscale*b  //bdbdphib
      *2.*constants::PI*p; //pdphip
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}

/////////////////////////////////////////////
//////// b dependent hadron <pt> ///////////
////////////////////////////////////////////

static int HadronsAvPtNum(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define havptnumk xx[0]
#define havptnumphik xx[1]
#define havptnumR xx[2]
#define havptnumphiR xx[3]
#define havptnumb xx[4]
#define havptnumphib xx[5]
#define havptnumphi xx[6]
#define havptnump xx[7]

  double z = xx[8];
  double mh = 0.2; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 2./constants::hbarc;
  double bscale = 12./constants::hbarc;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double lambda = static_cast<params*>(userdata)->lambda;
  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;
  double Y = static_cast<params*>(userdata)->Y;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double Bp = static_cast<params*>(userdata)->Bp;
  double A = static_cast<params*>(userdata)->A;
 
  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double phi = 2.*constants::PI*havptnumphi;
  double phik = 2.*constants::PI*havptnumphik;
  double phiR = 2.*constants::PI*havptnumphiR;
  double phib = 2.*constants::PI*havptnumphib;
  double k = havptnumk*kscale;
  double p = havptnump*pscale+lambda;
  double R = havptnumR*Rscale;
  double b = havptnumb*bscale;

  double eta = Y;
  double pg = p/z;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  // //double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  // double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;

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
    
    double Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    double QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    double TA = returnTA(sqrt(max(R*R + b*b - 2.*R*b*cos((phiR - phib)),0.)),TAclass);
    
    f = Dh*J*
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.)//!!
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      *Phip(k, R, Qsp, sizeFactor, mv, BK,xp,bdep,bdep_p, Bp, alphas)*factorxp
      *Phit(sqrt((pg)*(pg) + k*k - 2.*pg*k*cos((phi - phik))), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)
      *factorxA
      *2.*constants::PI*kscale*k  //kdkdphik
      *2.*constants::PI*Rscale*R  //RdRdphiR
      *2.*constants::PI*bscale*b  //bdbdphib
      *2.*constants::PI*pscale*p*p; //p^2dpdphip
    //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  }
  return 0;
}


/////////////////////////////////////////////////////
///// Fluctuating b hadrons cross section ///////////
////////////////////////////////////////////////////

static int HadronsFluc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

 #define z xx[6] 

  double mh = 0.2; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc;
  double bscale = 24./constants::hbarc;
  double Rx = fgRx*Rscale-Rscale/2.;
  double Ry = fgRy*Rscale-Rscale/2.;
  double k = fgk*kscale;

  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;

  double lambda = static_cast<params*>(userdata)->lambda;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double Y = static_cast<params*>(userdata)->Y;
  double eta = Y;
  double Qsp=0.;
  double QsA=0.;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double A = static_cast<params*>(userdata)->A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double R = sqrt(Rx*Rx+Ry*Ry);
  double b = sqrt(bx*bx+by*by);

  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
  
  double p = (fgp*pscale+lambda);
  
  double pg = p/z;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  // //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  // double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;
  
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
    
    Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
    f = Dh*J*
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2*constants::PI*constants::PI),3.)
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      *PhipFluc(k, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp
      *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos((fgphi - fgphik)*2.*constants::PI)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)*factorxA
      *2.*constants::PI*k*kscale  //kdkdphik
      *Rscale*Rscale  //dRxdRy
      //    *bscale*bscale  //bdxdby
      *2.*constants::PI*pscale*p; //pdpdphip
  }
  
  return 1;
}


static int HadronsFlucAvPtNum(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

 #define z xx[6] 

  double mh = 0.2; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc;
  double bscale = 24./constants::hbarc;
  double Rx = fgRx*Rscale-Rscale/2.;
  double Ry = fgRy*Rscale-Rscale/2.;
  double k = fgk*kscale;

  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;

  double lambda = static_cast<params*>(userdata)->lambda;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double Y = static_cast<params*>(userdata)->Y;
  double eta = Y;
  double Qsp=0.;
  double QsA=0.;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double A = static_cast<params*>(userdata)->A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double R = sqrt(Rx*Rx+Ry*Ry);
  double b = sqrt(bx*bx+by*by);
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
  
  double p = (fgp*pscale+lambda);
  
  double pg = p/z;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  // //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  // double Dh = kkp::KKPFragmentation(7, 1, z, p, gluon);
  
  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;
  
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
    
    Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
    f = Dh*J*
      alphas/constants::CF/(z*pg*z*pg+mIR*mIR)/pow((2*constants::PI*constants::PI),3.)
      //alphas/constants::CF/(z*z)/(pg*pg+mIR*mIR)/pow((2.*constants::PI*constants::PI),3.) //!!
      *PhipFluc(k, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp
      *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos((fgphi - fgphik)*2.*constants::PI)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)*factorxA
      *2.*constants::PI*k*kscale  //kdkdphik
      *Rscale*Rscale  //dRxdRy
      //    *bscale*bscale  //bdxdby
      *2.*constants::PI*pscale*p*p; //p^2dpdphip
  }
  
  return 1;
}

static int HadronsFlucNoPT(const int *ndim, const cubareal xx[],
                           const int *ncomp, cubareal ff[], void *userdata) {
  
#define fgk xx[0]
#define fgphik xx[1]
#define fgRx xx[2]
#define fgRy xx[3]
#define fgphi xx[4]
#define zfnpt xx[5] 

  double mh = 0.2; //GeV

  double kscale = 30.;
  double pscale = 30.;
  double Rscale = 4./constants::hbarc;
  double bscale = 24./constants::hbarc;
  double Rx = fgRx*Rscale-Rscale/2.;
  double Ry = fgRy*Rscale-Rscale/2.;
  double k = fgk*kscale;

  double bx=static_cast<params*>(userdata)->bx/constants::hbarc;
  double by=static_cast<params*>(userdata)->by/constants::hbarc;

  double lambda = static_cast<params*>(userdata)->lambda;
  double sizeFactor = static_cast<params*>(userdata)->protonSizeFactor;
  double Y = static_cast<params*>(userdata)->Y;
  double eta = Y;
  double Qsp=0.;
  double QsA=0.;
  double mIR = static_cast<params*>(userdata)->mIR;

  double alphas = static_cast<params*>(userdata)->alphas;
  double roots = static_cast<params*>(userdata)->roots;
  double sigma02 = static_cast<params*>(userdata)->sigma02;
  double rt2 = static_cast<params*>(userdata)->rt2;
  double bdep_p = static_cast<params*>(userdata)->bdep_p;
  double bindep_A = static_cast<params*>(userdata)->bindep_A;
  double bdep_A = static_cast<params*>(userdata)->bdep_A;
  double bdep_fluc_p = static_cast<params*>(userdata)->bdep_fluc_p;
  double bdep_fluc_A = static_cast<params*>(userdata)->bdep_fluc_A;
  double A = static_cast<params*>(userdata)->A;
  int Ybin = static_cast<params*>(userdata)->Ybin;

  int BK = static_cast<params*>(userdata)->BK;
  int bdep = static_cast<params*>(userdata)->bdep;
  int useFluc = static_cast<params*>(userdata)->useFluc;

  TAInt *TAclass = static_cast<params*>(userdata)->TAclass;
  Glauber *glauberClass = static_cast<params*>(userdata)->glauberClass;
  MV *mv = static_cast<params*>(userdata)->mv;

  double R = sqrt(Rx*Rx+Ry*Ry);
  double b = sqrt(bx*bx+by*by);
  double TA = returnTA2D(Rx-bx,Ry-by,glauberClass, Ybin);
  double Tp = returnTp2D(Rx,Ry,glauberClass, Ybin);
  
  double p = static_cast<params*>(userdata)->PT;
  
  double pg = p/zfnpt;
  
  //  double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  double J = p*cosh(eta)/sqrt(p*p*cosh(eta)*cosh(eta)+mh*mh);
  //  double Dh = 6.05*pow(z,-0.714)*pow(1.-z,2.92); //KKP NLO 
  double Dh = kkp::KKPFragmentation(7, 1, zfnpt, p, gluon);
  
  double yg;
  if (A==208)
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta)))) + 0.465;
  else
    yg = 0.5*log((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))+p*sinh(eta))
                 /((sqrt(mh*mh+p*p*cosh(eta)*cosh(eta))-p*sinh(eta))));

  // double J = pg*cosh(eta)/sqrt(pg*pg*cosh(eta)*cosh(eta)+mh*mh);
  // //  double Dh = 6.05*pow(zfnpt,-0.714)*pow(1.-zfnpt,2.92); //KKP NLO 
  // double Dh = kkp::KKPFragmentation(7, 1, zfnpt, p, gluon);
  
  // double yg;
  // if (A==208)
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta)))) + 0.465;
  // else
  //   yg = 0.5*log((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))+pg*sinh(eta))
  //                /((sqrt(mh*mh+pg*pg*cosh(eta)*cosh(eta))-pg*sinh(eta))));
  
  double xp = pg*exp(yg)/roots;
  double xA = pg*exp(-yg)/roots;
  
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
    
    Qsp = constants::prefactorp*pow(constants::x0/xp,constants::lambdaSpeedp/2.);
    QsA = constants::prefactorA*pow(constants::x0/xA,constants::lambdaSpeedA/2.);
    
    // Below use Phip(..,Tp,..) when using quarks in the proton, otherwise use Phip(..,R,..) 
    f = Dh*J* 
      //      alphas/constants::CF/(zfnpt*pg*zfnpt*pg+mIR*mIR)/pow((2*constants::PI*constants::PI),3.)
      alphas/constants::CF/(zfnpt*zfnpt)/(pg*pg+mIR*mIR)/pow((2*constants::PI*constants::PI),3.)
      *PhipFluc(k, Tp, Qsp, sizeFactor, mv, BK, xp, bdep_fluc_p, alphas)*factorxp
      *Phit(sqrt(pg*pg + k*k - 2.*pg*k*cos((fgphi - fgphik)*2.*constants::PI)), TA, QsA, mv, BK, xA,bdep,useFluc, bindep_A, bdep_A, bdep_fluc_A, alphas)*factorxA
      *2.*constants::PI*k*kscale  //kdkdphik
      *Rscale*Rscale  //dRxdRy
      //    *bscale*bscale  //bdxdby
      *2.*constants::PI*p; //2pi p
  }
  
  return 1;
}


// Main program
int main(int argc, char *argv[]) {
 // MPI things
  int rank=0;
  int size=1;

  // Option defaults
  int readTable = 1;
  int dopt = 0;
  int useFluc = 0;
  int Nevents = 1;
  int NRQCD = 1;
  int BK = 1;
  int BKMVe = 0;
  int bdep = 0;
  double Yg = 0.;
  double YJPsi1 = 2.73;
  double YJPsi2 = -3.71;
  int xsec = 0;
  double alphas = 0.3;
  string Target="Pb";
  double bmax = 10.;
  double Bp = 4.; // GeV^-2
  double Bq = 0.3;
  double sigma02In = 9.; // mb
  double RA = 5.5; // fm
  double width = 0.5; // Qs fluctutation sigma
  const double oomph = 2.21; //for Pb
  double mIR=0.2;
  double roots = 8160.;
  double Bqwidth = 0.;
  int Nq = 3;
  int BqYdep = 0;
  int flucNq = 0;

  // initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
  MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes

  params data;

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
        std::cerr << "--BK option requires one argument, 0 for MV with Qs(x), 1 for (at the moment fake) BK, or 2 for an interpolation of MV and BK (option 2 for testing in readTable mode only)." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--BKMVe") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        BKMVe = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--BKMVe option requires one argument, 0 for MV or 1 for MVe initial condition (only used if BK=1)" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--bdep") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        bdep = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--bdep option requires one argument, 0 for b-independent or 1 for b-dependent." << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Yg") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        Yg = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Yg option requires one argument, gluon rapidity" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--YJPsi1") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        YJPsi1 = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--YJPsi1 option requires one argument, first JPsi rapidity" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--YJPsi2") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        YJPsi2 = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--YJPsi2 option requires one argument, second JPsi rapidity" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--xsec") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        xsec = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--xsec option requires one argument, 0 or 1, with 1 meaning that we run an event to later determine the inelastic cross section" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--sigma02") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        sigma02In = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--sigma02 is sigma_0/2 in [mb], this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--RA") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        RA = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--RA is the radius of the target nucleus in fm; this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Bp") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        Bp = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Bp is the size of the proton in GeV^(-2); this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Bq") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        Bq = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Bq is the size of a hotspot in GeV^(-2); this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--alphas") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        alphas = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--alphas is the strong coupling constant; this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Qswidth") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        width = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Qswidth is the width of the Qs fluctuation distribution; this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--dopt") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        dopt = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--dopt selects if the pT spectrum should be calculated (1) or not (0) in the fluctuating case; this option requires one argument, 0 or 1" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Target") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        Target = argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Target selects the target nucleus, default is 'Pb'" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--bmax") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        bmax = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--bmax maximal impact parameter in [fm]; this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--mIR") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        mIR = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--mIR is IR regulator mass in [GeV] for the gluons and hadrons; this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--roots") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        roots = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--roots is the center of mass energy in [GeV]; this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Bqwidth") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        Bqwidth = atof(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Bqwidth is the width of the log-normal fluctuations of Bq; this option requires one argument, a double >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--Nq") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        Nq = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--Nq is the number of hot spots per nucleon; this option requires one argument, an integer >0" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--BqYdep") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        BqYdep = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--BqYdep decides if Bq should depend on rapidity; this option requires one argument, 0 or 1" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "--flucNq") {
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        flucNq = atoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--flucNq decides if Nq should fluctuate event to event (it is the same for all nucleons); this option requires one argument, 0 or 1" << std::endl;
        return 1;
      }  
    }
    else if (std::string(argv[i]) == "-?") {
      cout << "Options are:\n" << "--readTable [0 or 1], --fluctuations [0 or 1], --Nevents [# of events], --NRQCD [0 or 1], --BK [0 or 1], --BKMVe [0 or 1], --bdep [0 or 1], --Yg [value of gluon rapidity], --YJPsi1 [value of first JPsi rapidity], --YJPsi2 [value of second JPsi rapidity], --xsec [0 or 1] (computes inelastic cross section when set to 1, --sigma02 [>0. mb], --RA [>0. fm], --Bp [>0. GeV^-2], --Bq [>0. GeV^-2], --alphas [>0.], --Qswidth [>0.], --dopt [0 or 1], -- Target [Pb, p, ...], --bmax [>=0. fm], --mIR [>0. GeV], --roots [>0 GeV], --Bqwidth [>0.], --Nq [>0], --BqYdep [0 or 1])" << endl;
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i++;
        exit(0);
      } else {
        exit(0);
      }  
    }
  }



  Parameters *Glauber_param;
  Glauber_param = new Parameters();
  Glauber_param->setParameters(Target);
  
  
  long int seed = time(NULL)+rank*100000;
  //long int seed = 1;
    
  Random *random;
  random = new Random();
  random->init_genrand64(seed);
  
  Glauber *glauber;
  glauber = new Glauber(Glauber_param, width);
  glauber->init(random);
  double A = glauber->getA();

  data.roots = roots;
  data.mIR=mIR;
  data.alphas=alphas;
  data.Bp=Bp;
  data.Bq=Bq;
  data.sigma02 = sigma02In/10./constants::hbarc/constants::hbarc; 
  data.rt2 =  pow(RA,2.)/constants::hbarc/constants::hbarc;//(1.69/constants::hbarc/constants::hbarc)*pow(A,2./3.);// assuming R_A = R_0 A^(1/3), with R_0=1.3 fm (R_0^2 = 1.69 fm^2)
  if(useFluc == 0){
    data.bdep_p = data.sigma02/2./constants::PI/data.Bp;
  }
  else{
    data.bdep_p = data.sigma02/2./constants::PI/(data.Bp+data.Bq);
  }

  data.bdep_A = oomph*data.bdep_p;// this is for Pb
  //  const double bindep_A = 1./(oomph*(1.69/constants::hbarc/constants::hbarc)/pow(A,1./3.)/2./Bp)*bdep_A; // assuming R_A = R_0 A^(1/3), with R_0=1.3 fm (R_0^2 = 1.69 fm^2)
  data.bindep_A = A*data.sigma02/constants::PI/data.rt2; 
  data.bdep_fluc_p = data.bdep_p;
  data.bdep_fluc_A = data.bdep_p; // same as for the proton

  int h5Flag = 0;
  pretty_ostream messenger;
  
  
  if (rank==0){
    cout << "Options: \n read dipole from file yes(1)/no(0)= " << readTable 
         << "\n fluctuations on(1)/off(0) = " << useFluc 
         << "\n Number of events = " << Nevents 
         << "\n NRQCD on(1)/off(0) = " << NRQCD 
         << "\n BK on(1)/off(0) = " << BK
         << "\n MVe initial condition = " << BKMVe 
         << "\n b dependence on(1)/off(0) = " << bdep
         << "\n Y_g = " << Yg
         << "\n YJPsi1 = " << YJPsi1
         << "\n YJPsi2 = " << YJPsi2 
         << "\n cross section mode on(1)/off(0) = " << xsec 
         << "\n print pT spectrum (fluc. case) = " << dopt
         << endl;
    cout << " Target = " << Target << endl;
    cout << " A = " << A << endl;
    cout << " roots = " << roots << " GeV" << endl;
    cout << " mIR = " << mIR << " GeV" << endl;
    cout << " bmax = " << bmax << " fm" << endl;
    cout << " alphas = " << alphas << endl;
    cout << " sigma02 = " << sigma02In << " mb" << endl;
    cout << " sigma02 = " << data.sigma02 << " GeV^(-2)" << endl;
    cout << " RA = " << RA << " fm" << endl;
    cout << " rt2 (from RA) = " << data.rt2 << " GeV^(-2)" << endl;
    cout << " Nq = " << Nq << endl;
    cout << " Bp = " << data.Bp << " GeV^(-2)" << endl;
    cout << " Bq = " << data.Bq << " GeV^(-2)" << endl;
    cout << " Qswidth = " << width << endl;
    cout << " Bqwidth = " << Bqwidth << endl;
    cout << " BqYdep = "  << BqYdep << endl;
    cout << " flucNq = "  << flucNq << endl;
    cout << " bdep_p = " << data.bdep_p << endl;
    cout << " bdep_A = " << data.bdep_A << endl;
    cout << " bindep_A = " << data.bindep_A << endl;
    
    messenger.flush("info");
  }

  

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
      mv->computePhipBK(BKMVe);
      mv->writeTableBK(BKMVe);
    }
  }
  else if (readTable == 1){
    if(BK==0){
      mv->readTable();
    }
    else if(BK==1){
      mv->readTableBK(BKMVe);
    }
    else if(BK==2){
      mv->readTable();
      mv->readTableBK(BKMVe);
    }
  }
  else{
    cerr << "Unknown option readTable = " << readTable << ". Options are 0 or 1." << endl;
    exit(0);
  }


 //  //test Phit
 
 //  stringstream strfilenamet;
 //  strfilenamet << "Phit.dat";
 //  string filenamet;
 //  filenamet = strfilenamet.str();
 //  fstream foutt(filenamet.c_str(), ios::out);
 //  double QsIn = constants::prefactor*pow(constants::x0/0.00001,constants::lambdaSpeedp/2.);
 //  //  cout << "Qs=" << QsIn << endl;
 //  for (int ik=0; ik<1000; ik++){
 //   double k = ik*0.01;
 //   //    double StF(double k, double TA, double Qs, MV *mv, int BK, double x, int bdep, int useFluc){
 //   foutt << k << " " << Phit(k, 1., QsIn, mv, BK, 0.00001, 1, 0)  << endl;
 // }


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
 
  int comp, nregions, fail;
  long long int neval;

  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  double gresult = 0.;
  double gerror = 0.;

  double hresult = 0.;
  double herror = 0.;
  double hresultpt = 0.;
  
  double hmeanPt = 0;
  double hmeanPterror_num = 0;
  double hmeanPtresult_num = 0;
  double hmeanPtresult_den = 0;
  double hmeanPterror_den = 0;

  double JPsi2result = 0.;
  double JPsi2error = 0.;

  double JPsi2result_cs = 0.;
  double JPsi2error_cs = 0.;

  double JPsi2result_co = 0.;
  double JPsi2error_co = 0.;

  double JPsi2result_num = 0.;
  double JPsi2error_num = 0.;

  double JPsi2result_den = 0.;
  double JPsi2error_den = 0.;

  double JPsi2result_num_2 = 0.;
  double JPsi2error_num_2 = 0.;

  double JPsi2result_den_2 = 0.;
  double JPsi2error_den_2 = 0.;

  double JPsi2result2 = 0.;
  double JPsi2error2 = 0.;

  double meanPtresult_num = 0;
  double meanPterror_num = 0;

  double meanPtresult_den = 0;
  double meanPterror_den = 0;

  double Jpsi2meanPtresult = 0;
  double Jpsi2meanPtresult_2 = 0;

  data.useFluc = useFluc;
  data.bdep = bdep;
  data.BK = BK;

  data.PT = 0.; // dummy for now
  data.pe = 0.; // dummy for now
  data.k = 0.;  // dummy for now
  data.m = 0.;
  data.lambda = 0.05; // Infrared cutoff on p integral in GeV (50 MeV according to https://arxiv.org/pdf/1812.01312.pdf) - use 0 when using mIR as cutoff
  data.mv = mv; // MV class
  data.TAclass = TAclass; // TA class
  data.glauberClass = glauber; // Glauber class
  data.protonSizeFactor = 1.;
  data.A = A;

  data.Y = Yg;  

  if(useFluc == 0){
    cout << "For b integrated results obtained in this mode (no fluctuations) all results are cross sections, that need to be divided by the total inelastic cross section (in p+Pb) to get particle numbers." << endl; 
    if(bdep == 0){
      cout << "b-independent results"  << endl;
 
      // Code to compute the pt spectrum
      int npoints = 20;
      
      for (int nip=0; nip<npoints; nip++){ 
        //       data.PT= pow(10,ptmin + nip*step);   
        data.PT= 0.1 + (nip*0.5);
        NDIM = 4;
        llVegas(NDIM, NCOMP, HadronsNoBNoPt, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hresult = (double)integral[0];
        NDIM = 7;
        llVegas(NDIM, NCOMP, GluonsNoBNoPt, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        gresult = (double)integral[0];
        NDIM = 11;
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDNoBNoPt, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        JPsi2result = (double)integral[0];
        
        cout << data.PT << " " << hresult << " " << gresult << " " << JPsi2result << endl;
      }
      
      cout << " -- -- -- -- -- -- -- -- --" << endl;



       //for (int nip=0; nip<=npoints; nip++){ 
        //   data.PT= pow(10,ptmin + nip*step);
         //  NDIM = 4;
       //  llVegas(NDIM, NCOMP, HadronsNoBNoPt, &data, NVEC,
       //         EPSREL, EPSABS, VERBOSE, SEED,
        //        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        //        GRIDNO, NULL, NULL,
        //        &neval, &fail, integral, error, prob);
      //}
      
      for(int i=0; i<9;i++){
        data.Y = -4.+i;
        NDIM = 4;
        llVegas(NDIM, NCOMP, GluonsNoB, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        gresult = (double)integral[0];
        gerror = (double)error[0];

        NDIM = 5;
        llVegas(NDIM, NCOMP, HadronsNoB, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hresult = (double)integral[0];
        herror = (double)error[0];
        
        data.lambda = 0.15; // Infrared cutoff on p integral in GeV (50 MeV according to https://arxiv.org/pdf/1812.01312.pdf)
        NDIM = 5;
        llVegas(NDIM, NCOMP, HadronsNoBAvPtNum, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hmeanPtresult_num = (double)integral[0];
        hmeanPterror_num = (double)error[0];

        llVegas(NDIM, NCOMP, HadronsNoB, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hmeanPtresult_den = (double)integral[0];
        hmeanPterror_den = (double)error[0];
        
        hmeanPt = hmeanPtresult_num/hmeanPtresult_den;
        data.lambda = 0.05; // Infrared cutoff on p integral in GeV (50 MeV according to https://arxiv.org/pdf/1812.01312.pdf)
       
        if(NRQCD==1){
          //    cout << "Using NRQCD"  << endl;
          NDIM = 8;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsNoB, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          
          JPsi2result_cs = (double)integral[0];
          JPsi2error_cs = (double)error[0];    
          
          NDIM = 6;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoNoB, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          JPsi2result_co= (double)integral[0];
          JPsi2error_co = (double)error[0];
          
          JPsi2result = JPsi2result_cs + JPsi2result_co;
          JPsi2error = JPsi2error_cs + JPsi2error_co;

          NDIM = 8;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDNoBAvPtNum, &data, NVEC,
                 EPSREL, EPSABS, VERBOSE, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, NULL, NULL,
                 &neval, &fail, integral, error, prob);

         meanPtresult_num= (double)integral[0];
         meanPterror_num = (double)error[0];

          NDIM = 8;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDNoBAvPtDen, &data, NVEC,
                 EPSREL, EPSABS, VERBOSE, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, NULL, NULL,
                 &neval, &fail, integral, error, prob);

         meanPtresult_den = (double)integral[0];
         meanPterror_den = (double)error[0];

         Jpsi2meanPtresult = meanPtresult_num/meanPtresult_den;

        }
        else{
          //cout << "Using ICEM"  << endl;    
          NDIM = 8;
          llVegas(NDIM, NCOMP, JPsiIntegrandAllNoB, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          JPsi2result = (double)integral[0];
          JPsi2error = (double)error[0];  
        }
        cout << data.Y << " " << setprecision(10) << gresult << " " << gerror << " " 
             << hresult << " " << herror << " " << hmeanPt << " " << JPsi2result
             << " " << JPsi2error << " " << Jpsi2meanPtresult << endl;
      }
    }
    else{
      cout << "b-dependent results"  << endl;

      // Code to compute the pt spectrum
      int npoints = 20;
      
      for (int nip=0; nip<npoints; nip++){ 
        //       data.PT= pow(10,ptmin + nip*step);   
        data.PT= 0.1 + (nip*0.5);
        NDIM = 8;
        llVegas(NDIM, NCOMP, HadronsNoPt, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hresult = (double)integral[0];
        NDIM = 7;
        llVegas(NDIM, NCOMP, GluonsNoPt, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        gresult = (double)integral[0];
        NDIM = 11;
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDNoPT, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        JPsi2result = (double)integral[0];
        
        cout << data.PT << " " << hresult << " " << gresult << " " << JPsi2result << endl;
      }
      
      cout << " -- -- -- -- -- -- -- -- --" << endl;

      for(int i=0; i<9;i++){
        data.Y = -4.+i;
        
        NDIM = 8;
        llVegas(NDIM, NCOMP, Gluons, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        gresult = (double)integral[0];
        gerror = (double)error[0];

        NDIM = 9;
        llVegas(NDIM, NCOMP, Hadrons, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hresult = (double)integral[0];
        herror = (double)error[0];

        data.lambda = 0.15; // Infrared cutoff on p integral in GeV (50 MeV according to https://arxiv.org/pdf/1812.01312.pdf)
        NDIM = 9;
        llVegas(NDIM, NCOMP, Hadrons, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hmeanPtresult_den = (double)integral[0];
        hmeanPterror_den = (double)error[0];
        NDIM = 9;
        llVegas(NDIM, NCOMP, HadronsAvPtNum, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        hmeanPtresult_num = (double)integral[0];
        hmeanPterror_num = (double)error[0];

        hmeanPt = hmeanPtresult_num/hmeanPtresult_den;
        data.lambda = 0.05; 
         
        if(NRQCD==1){
          //cout << "Using NRQCD"  << endl;
          NDIM = 12;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCs, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          
          JPsi2result_cs = (double)integral[0];
          JPsi2error_cs = (double)error[0];    
          
          NDIM = 10;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCo, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          JPsi2result_co= (double)integral[0];
          JPsi2error_co = (double)error[0];
          
          JPsi2result = JPsi2result_cs + JPsi2result_co;
          JPsi2error = JPsi2error_cs + JPsi2error_co;

          NDIM = 12;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtNum, &data, NVEC,
                 EPSREL, EPSABS, VERBOSE, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, NULL, NULL,
                 &neval, &fail, integral, error, prob);

         meanPtresult_num= (double)integral[0];
         meanPterror_num = (double)error[0];

          NDIM = 12;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtDen, &data, NVEC,
                 EPSREL, EPSABS, VERBOSE, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, NULL, NULL,
                 &neval, &fail, integral, error, prob);

         meanPtresult_den = (double)integral[0];
         meanPterror_den = (double)error[0];

         Jpsi2meanPtresult = meanPtresult_num/meanPtresult_den;
 
          //          cout << data.Y << " " << gresult << " " << hresult << " " << JPsi2result_cs+JPsi2result_co << endl;
        }
        else{
          //cout << "Using ICEM"  << endl;    
          NDIM = 12;
          llVegas(NDIM, NCOMP, JPsiIntegrandAll, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          JPsi2result = (double)integral[0];
          JPsi2error = (double)error[0];  
          //cout << data.Y << " " << gresult << " " << hresult << " " << JPsi2result << endl;
        }
        cout << data.Y << " " << setprecision(10) << gresult << " " << gerror << " " 
             << hresult << " " << herror << " " << hmeanPt << "  " << JPsi2result
             << " " << JPsi2error << " " << Jpsi2meanPtresult << endl;
  
      }
    }
      
  }
  
  else{
    if (rank==0)
      cout << "Fluctuating results"  << endl; 
    int ni=0;
    default_random_engine generator;    
    generator.seed(seed);    
    //double BqGauss;
    //double BqGaussSum =0.;
    //double BqGaussSumSq = 0.;
    while (ni<Nevents){
      
      // fluctuation of the hot-spot size (moved to Glauber.cpp to fluctuate each hot spot differently
      // for (int i=0;i<10000;i++){
      //   BqGaussSum += BqGauss;
      //   BqGaussSumSq += BqGauss*BqGauss;
      // }
      // cout << "E[Gauss]=" << BqGaussSum/10000. << "sigma[Gauss]=" << sqrt(BqGaussSumSq/10000.-BqGaussSum*BqGaussSum/10000./10000.) << endl;
      
      // Run Vegas integration with fluctuations
      // Make a new target
      
      if( flucNq == 1 ){
        poisson_distribution<int> d(2.8214393721220787); // this mean makes a mean of 3 for the zero truncated Poisson distribution
        Nq = d(generator);
        while(Nq<1){
          Nq = d(generator);
        }
        //cout << "Nq = " << Nq << endl;
      }


      glauber->makeNuclei(random, data.Bp, data.Bq, Bqwidth, Nq, Yg, YJPsi1, YJPsi2, BqYdep);
      
      // Sample b
      double bmin = 0.;
      
      double xb =
        random->genrand64_real1(); // uniformly distributed random variable
      double b = sqrt((bmax * bmax - bmin * bmin) * xb + bmin * bmin);
      double phib = 2.*constants::PI*random->genrand64_real1();
      
      data.bx = b*cos(phib);
      data.by = b*sin(phib);
      data.Y = Yg;
      data.Ybin = 0;
      if(BqYdep){
        data.bdep_fluc_p = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
        data.bdep_fluc_A = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((-data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
      }
      
      // do hadrons
      NDIM = 7;
      llVegas(NDIM, NCOMP, HadronsFluc, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);
      
      // Print the result
      hresult = (double)integral[0];
      herror = (double)error[0];
      //printf("Hadrons (fluc): %.8f +- %.8f\t\n", hresult, herror);
      
      if(xsec == 1){
        if(hresult<0.5){
          //  cout << "Hadron number < 0.5, skipping event" << endl;
          exit(0);
        }
        else{
          stringstream strfilenameh;
          strfilenameh << "output_h_" << rank << ".dat";
          string filenameh;
          filenameh = strfilenameh.str();
          fstream fouth(filenameh.c_str(), ios::app);
          
          fouth << std::scientific << setprecision(5) << hresult << " " << herror << " " << 0.
                << " " << 0. << " " << 0. << " " << 0. << endl;
          fouth.close();
          exit(0);
        }
      }

      // the data (at least in pp) requires at least one charged particle track in |eta|<1. (https://arxiv.org/pdf/1202.2816.pdf)
      if(hresult<0.5){
        //        cout << "Hadron number < 0.5, skipping event" << endl;
        continue;
      }

      if (dopt==1){
        stringstream strfilenamehpt;
        strfilenamehpt << "output_hpt_Yg" << Yg << "_" << rank << "_" << ni << ".dat";
        //strfilenameh << "output_h_" << rank << ".dat";
        string filenamehpt;
        filenamehpt = strfilenamehpt.str();
        fstream fouthpt(filenamehpt.c_str(), ios::out);
        
        // Code to compute the pt spectrum
        int npoints = 20;
        
        for (int nip=0; nip<npoints; nip++){ 
          //       data.PT= pow(10,ptmin + nip*step);   
          data.PT= 0.1 + (nip*0.5);
          data.Y = Yg;
          data.Ybin = 0;
          NDIM = 6;
          llVegas(NDIM, NCOMP, HadronsFlucNoPT, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          hresultpt = (double)integral[0];
          data.Y = YJPsi1; //forward
          data.Ybin = 1;
          NDIM = 9;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDFlucNoPT, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          JPsi2result = (double)integral[0];
          // llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsFlucNoPT, &data, NVEC,
          //         EPSREL, EPSABS, VERBOSE, SEED,
          //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          //         GRIDNO, NULL, NULL,
          //         &neval, &fail, integral, error, prob);
          // JPsi2result_cs = (double)integral[0];
          // NDIM = 7;
          // llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoFlucNoPT, &data, NVEC,
          //         EPSREL, EPSABS, VERBOSE, SEED,
          //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          //         GRIDNO, NULL, NULL,
          //         &neval, &fail, integral, error, prob);
          // JPsi2result_co = (double)integral[0];
          
          cout << data.PT << " " << hresultpt << " " << JPsi2result << endl;
          fouthpt << std::scientific << setprecision(5) << data.PT << " " << hresultpt << " " << JPsi2result << endl;
        }
        fouthpt.close();
      }
      // done pT spectra

      data.Y = Yg;
      data.Ybin = 0;
      NDIM = 6;
      llVegas(NDIM, NCOMP, GluonsFluc, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);
      gresult = (double)integral[0];
      gerror = (double)error[0];
      //      printf("Midrapidity gluon (fluc): %.8f +- %.8f\t\n", gresult, gerror);        

      data.lambda = 0.15; // Infrared cutoff on p integral in GeV (50 MeV according to https://arxiv.org/pdf/1812.01312.pdf)
      NDIM = 7;
      llVegas(NDIM, NCOMP, HadronsFluc, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);
      hmeanPtresult_den = (double)integral[0];
      hmeanPterror_den = (double)error[0];
      NDIM = 7;
      llVegas(NDIM, NCOMP, HadronsFlucAvPtNum, &data, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, NULL, NULL,
              &neval, &fail, integral, error, prob);
      hmeanPtresult_num = (double)integral[0];
      hmeanPterror_num = (double)error[0];
      
      hmeanPt = hmeanPtresult_num/hmeanPtresult_den;
      data.lambda = 0.05; 
      
      if(NRQCD==1){
        //        cout << "Using NRQCD"  << endl;
        
        data.Y = YJPsi1; //forward
        data.Ybin = 1;
        if(BqYdep){
          data.bdep_fluc_p = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
          data.bdep_fluc_A = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((-data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
        }

        // NDIM = 10;
        // llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsFluc, &data, NVEC,
        //         EPSREL, EPSABS, VERBOSE, SEED,
        //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        //         GRIDNO, NULL, NULL,
        //         &neval, &fail, integral, error, prob);
        
        // JPsi2result_cs = (double)integral[0];
        // JPsi2error_cs = (double)error[0];    
        
        // NDIM = 8;
        // llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoFluc, &data, NVEC,
        //         EPSREL, EPSABS, VERBOSE, SEED,
        //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        //         GRIDNO, NULL, NULL,
        //         &neval, &fail, integral, error, prob);
        // JPsi2result_co= (double)integral[0];
        // JPsi2error_co = (double)error[0];
        
        // JPsi2result = JPsi2result_co + JPsi2result_cs;
        // JPsi2error = JPsi2error_co + JPsi2error_cs;

        //mean pT
        NDIM = 10;
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtNumFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        
        JPsi2result_num = (double)integral[0];
        JPsi2error_num = (double)error[0];    
        
        llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtDenFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        
        JPsi2result_den = (double)integral[0];
        JPsi2error_den = (double)error[0];    
  
        Jpsi2meanPtresult = JPsi2result_num/JPsi2result_den;
        JPsi2result = JPsi2result_den;
        JPsi2error = JPsi2error_den;
          
        // NDIM = 12;
        // llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtNum, &data, NVEC,
        //          EPSREL, EPSABS, VERBOSE, SEED,
        //          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        //          GRIDNO, NULL, NULL,
        //          &neval, &fail, integral, error, prob);

        //  meanPtresult_num= (double)integral[0];
        //  meanPterror_num = (double)error[0];

        //   NDIM = 12;
        //   llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtDen, &data, NVEC,
        //          EPSREL, EPSABS, VERBOSE, SEED,
        //          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        //          GRIDNO, NULL, NULL,
        //          &neval, &fail, integral, error, prob);

        //  meanPtresult_den = (double)integral[0];
        //  meanPterror_den = (double)error[0];

        //  Jpsi2meanPtresult = meanPtresult_num/meanPtresult_den;

        if ( YJPsi1 == YJPsi2 ){
          JPsi2result2 = JPsi2result;
          JPsi2error2 = JPsi2error;
          Jpsi2meanPtresult_2 = Jpsi2meanPtresult;
        }
        else{
          data.Y = YJPsi2; //backward
          data.Ybin = 2;
          if(BqYdep){
            data.bdep_fluc_p = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
            data.bdep_fluc_A = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((-data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
          }
          
          // NDIM = 10;
          // llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCsFluc, &data, NVEC,
          //         EPSREL, EPSABS, VERBOSE, SEED,
          //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          //         GRIDNO, NULL, NULL,
          //         &neval, &fail, integral, error, prob);
          
          // JPsi2result_cs = (double)integral[0];
          // JPsi2error_cs = (double)error[0];    
          
          // NDIM = 8;
          // llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDCoFluc, &data, NVEC,
          //         EPSREL, EPSABS, VERBOSE, SEED,
          //         MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          //         GRIDNO, NULL, NULL,
          //         &neval, &fail, integral, error, prob);
          // JPsi2result_co= (double)integral[0];
          // JPsi2error_co = (double)error[0];
          
          //mean pT
          NDIM = 10;
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtNumFluc, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          
          JPsi2result_num_2 = (double)integral[0];
          JPsi2error_num_2 = (double)error[0];   
          
          llVegas(NDIM, NCOMP, JPsiIntegrandNRQCDAvPtDenFluc, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          
          JPsi2result_den_2 = (double)integral[0];
          JPsi2error_den_2 = (double)error[0];    
          
          Jpsi2meanPtresult_2 = JPsi2result_num_2/JPsi2result_den_2;
          
          JPsi2result2 = JPsi2result_den_2;
          JPsi2error2 = JPsi2error_den_2;
        }
        cout << setprecision(10) << hresult << " " << herror << " " << JPsi2result
             << " " << JPsi2error << " " << JPsi2result2 << " " << JPsi2error2
             << " " << hmeanPt << " " << Jpsi2meanPtresult << " " << Jpsi2meanPtresult_2 << endl;
      }
      else{
        //cout << "Using ICEM"  << endl;    

        data.Y = YJPsi1; //forward
        data.Ybin = 1;
        
        if(BqYdep){
          data.bdep_fluc_p = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
          data.bdep_fluc_A = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((-data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
        }

        //mean pT
        NDIM = 10;
        llVegas(NDIM, NCOMP, JPsiIntegrandICEMNumFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        JPsi2result_num = (double)integral[0];
        JPsi2error_num = (double)error[0];  

        llVegas(NDIM, NCOMP, JPsiIntegrandICEMDenFluc, &data, NVEC,
                EPSREL, EPSABS, VERBOSE, SEED,
                MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                GRIDNO, NULL, NULL,
                &neval, &fail, integral, error, prob);
        JPsi2result_den = (double)integral[0];
        JPsi2error_den = (double)error[0];  
  
        Jpsi2meanPtresult = JPsi2result_num/JPsi2result_den;
        JPsi2result = JPsi2result_den;
        JPsi2error = JPsi2error_den;
          
        if ( YJPsi1 == YJPsi2 ){
          JPsi2result2 = JPsi2result;
          JPsi2error2 = JPsi2error;
          Jpsi2meanPtresult_2 = Jpsi2meanPtresult;
        }
        else{
          data.Y = YJPsi2; //backward
          data.Ybin = 2;

          if(BqYdep){
            data.bdep_fluc_p = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
            data.bdep_fluc_A = data.sigma02/2./constants::PI/(data.Bp+data.Bq*(0.15 + 0.042*pow((-data.Y - 4.6),2.))); //make sure to use the same parametrization as in Glauber.cpp
          }
                    
          NDIM = 10;
          llVegas(NDIM, NCOMP, JPsiIntegrandICEMNumFluc, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);
          
          JPsi2result_num_2 = (double)integral[0];
          JPsi2error_num_2 = (double)error[0];   
          
          llVegas(NDIM, NCOMP, JPsiIntegrandICEMDenFluc, &data, NVEC,
                  EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                  GRIDNO, NULL, NULL,
                  &neval, &fail, integral, error, prob);

          JPsi2result_den_2 = (double)integral[0];
          JPsi2error_den_2 = (double)error[0];              
          
          Jpsi2meanPtresult_2 = JPsi2result_num_2/JPsi2result_den_2;
          
          JPsi2result2 = JPsi2result_den_2;
          JPsi2error2 = JPsi2error_den_2;
        }  

        cout << setprecision(10) << hresult << " " << herror << " " << JPsi2result
             << " " << JPsi2error << " " << JPsi2result2 << " " << JPsi2error2
             << " " << hmeanPt << " " << Jpsi2meanPtresult << " " << Jpsi2meanPtresult_2 << endl;

      }
      
      stringstream strfilename;
      strfilename << "output_g_Yg" << Yg << "_" << rank << ".dat";
      //strfilename << "output_g_" << rank << ".dat";
      string filename;
      filename = strfilename.str();
      fstream fout(filename.c_str(), ios::app);
      
      fout << std::scientific << setprecision(5) << gresult << " " << gerror << " " << JPsi2result
           << " " << JPsi2error << " " << JPsi2result2 << " " << JPsi2error2 << " "
           << sqrt(data.bx*data.bx+data.by*data.by)
           << endl;
      fout.close();

      stringstream strfilenameh;
      strfilenameh << "output_h_Yg" << Yg << "_" << rank << ".dat";
      //strfilenameh << "output_h_" << rank << ".dat";
      string filenameh;
      filenameh = strfilenameh.str();
      fstream fouth(filenameh.c_str(), ios::app);
      
      fouth << std::scientific << setprecision(5) << hresult << " " << herror << " " << JPsi2result
            << " " << JPsi2error << " " << JPsi2result2 << " " << JPsi2error2 << " " <<  hmeanPt << " " 
            << Jpsi2meanPtresult << " " 
            << Jpsi2meanPtresult_2 << endl;
      fouth.close();
     
      ni++;
  
    }
  }

  delete Glauber_param;
  delete random;
  delete mv;
  delete TAclass;
  delete glauber;

  MPI_Finalize();

  return 1;

}
