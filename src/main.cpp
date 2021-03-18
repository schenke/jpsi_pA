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
#include <vector>

#include <gsl/gsl_integration.h>
#include "cuba.h"
#include "pretty_ostream.h"
#include "hard.h"



#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

void display_logo();

// Constants to be used
namespace constants {
  const int Nc = 3;
  const double hbarc = 0.1973269804;
  const double CA = double(Nc);
  const double CF = (double(Nc)*double(Nc) - 1.)/(2.*double(Nc));
  const double alphas = 0.3;
  const double Bp = 0.4;
  const double mD = 1.864;
  const double mc = 1.275; //vary?
  const double mJPsi =3.096916;
}

// Parameters that need to be passed to the integrand
struct params {
  double pe;
  double k;
  double Qs;
  double lambda;
  double Y;
  double m;
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
  

// Nuclear density for lead
double rhoA(double z, void * params) {
  // do everything in GeV^some power
  double RA = 6.62/constants::hbarc;
  double d = 0.546/constants::hbarc;
  double R = *(double *) params;
  double f = 1/(1 + exp((sqrt(R*R + z*z) - RA)/d))/1.;
  return f;
}

double TA(double R){
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;
  double expected = -4.0;

  gsl_function F;
  F.function = &rhoA;
  F.params = &R;

  gsl_integration_qagi(&F, 1e-12, 1e-7, 1000, w, &result, &error);

  gsl_integration_workspace_free (w);
  
  return result/67.09678472225216694;// normalization for above parameters RA=6.62fm and d=0.546fm - adjust if parameters change;
}

// Unintegrated gluon distribution for the proton
double Phip(double k, double R, double Qs){
  return k*k*constants::CF/(pow(2.*M_PI,3))/constants::alphas
    *2.*constants::CF*exp(R*R/(2.*constants::Bp))*exp(-constants::CF
                                                     *exp(R*R/(2.*constants::Bp))*k*k/(constants::CA*Qs*Qs))/(constants::CA*Qs*Qs);
}

// Unintegrated gluon distribution for lead
double Phit(double k, double R, double Qs){
  return k*k*constants::CF/(pow(2.*M_PI,3))/constants::alphas*2.*constants::CF*exp(-constants::CF*k*k/(constants::CA*Qs*Qs*TA(R)))/(constants::CA*Qs*Qs*TA(R));
}

// FT of fundamental S of the target
double StF(double k, double R, double Qs){
  return 2.*exp(-(k*k/(Qs*Qs*TA(R))))/(Qs*Qs*TA(R));
} 

// Integrand for the first J/Psi integral
static int JPsiIntegrand1(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define qqR xx[0]
#define qqphiR xx[1]
#define qqb xx[2]
#define qqphib xx[3]
#define qqqtilde xx[4]
#define qqphi xx[5]
#define qqM xx[6]
#define qqPT xx[7]
#define qqphiPT xx[8]
#define f ff[0]

  double pscale = 100.;
  double Rscale = 3./constants::hbarc; //choose a small scale (proton Phip will cut off at large R)
  double bscale = 10./constants::hbarc; // bscale needs to be the same in all terms
  // Qs will be made rapidity dependent
  double Qs = static_cast<params*>(userdata)->Qs;
  double Y = static_cast<params*>(userdata)->Y;
  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + qqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = qqqtilde*qtildescale;
  double PT = qqPT*pscale; 
  double phiPT = qqphiPT*2.*M_PI;
  double R = qqR*Rscale;
  double b = qqb*Rscale;
  double phiR = qqphiR*2*M_PI;
  double phib = qqphib*2*M_PI;
  
  kinPair in;
  in.M = M;
  in.PT  = PT;
  in.phi = qqphi*2.*M_PI; // not the PT phi
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
  
  double pplusqx = px+qx;
  double pplusqy = py+qy;
  
  double pplusq = sqrt(pplusqx*pplusqx+pplusqy*pplusqy);
  double phi_pplusq = atan2(pplusqy,pplusqx);
  
  // since these 3 parts only appear in a sum, we should combine them to do only one
  // function call
  double H = Hard::qqqq(p, phip, q, phiq, pplusq, phi_pplusq, 0., 0., 0., 0., yp, yq, m)
    +Hard::qqg(p, phip, q, phiq, pplusq, phi_pplusq, 0., 0., 0., 0., yp, yq, m)
    +Hard::gg(p, phip, q, phiq, pplusq, phi_pplusq, 0., 0., 0., 0., yp, yq, m);
  
  // get Jacobian
  double betax = PT/sqrt(M*M+PT*PT);
  double gammax = 1./sqrt(1.-betax*betax);
  double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));

  f = constants::alphas*double(constants::Nc)*double(constants::Nc)
    /(2.*pow(2.*M_PI,8.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(pplusq, R, Qs)/(pplusq*pplusq)*H*J*(M/constants::mJPsi)*(M/constants::mJPsi)
    *R*Rscale*2.*M_PI
    *b*bscale*2.*M_PI
    *PT*pscale*2.*M_PI
    *(2.*constants::mD-constants::mJPsi)
    *qtildescale
    *2.*M_PI;
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // d2PT
  // dM
  // dqtilde
  // dphi

  // if (f < 0.)
  //   {
  //     cout << H << endl;
  //   }

  //remember factor 2 for p and q direction (not included yet)

  return 0;
}  
 
// Integrand for the second J/Psi integral
static int JPsiIntegrand2(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

  // others defined in Integrand1
#define qqk1 xx[9]
#define qqphik1 xx[10]

  double kscale = 100.;
  double pscale = 10.;
  double Rscale = 12./constants::hbarc;
  // Qs will be made rapidity dependent
  double Qs = static_cast<params*>(userdata)->Qs;
  double Y = static_cast<params*>(userdata)->Y;
  double m = constants::mc;//static_cast<params*>(userdata)->m;

  // scale the integration variables 
  double M = constants::mJPsi + qqM*(2.*constants::mD-constants::mJPsi); 
  double qtildescale = sqrt(M*M/4.-constants::mc*constants::mc);
  double qtilde = qqqtilde*qtildescale;
  double PT = qqPT*pscale; 
  double phiPT = qqphiPT*2.*M_PI;
  double R = qqR*Rscale;
  double b = qqb*Rscale;
  double k1 = qqk1*kscale;
  double phik1 = qqphik1*2.*M_PI;
  double phiR = qqphiR*2*M_PI;
  double phib = qqphib*2*M_PI;

  kinPair in;
  in.M = M;
  in.PT  = PT;
  in.phi = qqphi*2.*M_PI; // not the PT phi
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
  double k1x = k1*cos(phik1);
  double k1y = k1*sin(phik1);
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
  
  // since these 3 parts only appear in a sum, we should combine them to do only one
  // function call
  double H = Hard::qqqq(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, 0, 0, yp, yq, m)
    +Hard::qqg(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, 0, 0, yp, yq, m)
    +Hard::gg(p, phip, q, phiq, k1, phik1, pplusqminusk1, phi_pplusqminusk1, 0, 0, yp, yq, m);
  
  // get Jacobian
  double betax = PT/sqrt(M*M+PT*PT);
  double gammax = 1./sqrt(1.-betax*betax);
  double J = qtilde*gammax/(sqrt(p*p+m*m)*sqrt(q*q+m*m)*abs(sinh(yp-yq)));

  f = -constants::alphas*double(constants::Nc)*double(constants::Nc)
    /(2.*pow(2.*M_PI,9.)*(double(constants::Nc)*double(constants::Nc)-1.))
    *Phip(k1, R, Qs)/(k1*k1)*H*J*StF(pplusqminusk1,Rminusb,Qs)
    *R*Rscale*2.*M_PI
    *b*Rscale*2.*M_PI
    *PT*pscale*2.*M_PI
    *(2.*constants::mD-constants::mJPsi)
    *qtildescale
    *2.*M_PI
    *k1*kscale*2.*M_PI;
  // scaled momenta above (in PT)
  // last rows are scaling of integration measures:
  // d2R
  // d2b
  // d2PT
  // dM
  // dqtilde
  // dphi
  // d2k1
  //remember factor 2 for p and q direction (not included yet)

  return 0;
}  
  




// Integrand for integral over everything but |p|
static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define k xx[0]
#define phik xx[1]
#define R xx[2]
#define phiR xx[3]
#define b xx[4]
#define phib xx[5]
#define phi xx[6]

  double kscale = 100.;
  double Rscale = 12./constants::hbarc;
  double p = static_cast<params*>(userdata)->pe;
  double Qs = static_cast<params*>(userdata)->Qs;

  f = 2*M_PI*k*kscale*R*Rscale*b*Rscale*Phip(k*kscale, R*Rscale, Qs)*Phit(sqrt(p*p + k*k*kscale*kscale - 2.*p*k*kscale*cos((phi - phik)*2.*M_PI)), sqrt(max(R*Rscale*R*Rscale + b*b*Rscale*Rscale - 2.*R*b*Rscale*Rscale*cos((phiR - phib)*2.*M_PI),0.)), Qs)*2*M_PI*2*M_PI*kscale*Rscale*2*M_PI*Rscale*2*M_PI; // bscale = Rscale //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube)
  return 0;
}

// Integrand for the full 8D integral
static int FullIntegrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define p xx[7]

  double kscale = 100.;
  double pscale = 10.;
  double Rscale = 12./constants::hbarc;
  double Qs = static_cast<params*>(userdata)->Qs;
  double lambda = static_cast<params*>(userdata)->lambda;

  f = 2.*constants::alphas/constants::CF/(p*pscale+lambda)/(p*pscale+lambda)*2.*M_PI*k*kscale*R*Rscale*b*Rscale*Phip(k*kscale, R*Rscale, Qs)*Phit(sqrt((p*pscale+lambda)*(p*pscale+lambda) + k*k*kscale*kscale - 2.*(p*pscale+lambda)*k*kscale*cos((phi - phik)*2.*M_PI)), sqrt(max(R*Rscale*R*Rscale + b*b*Rscale*Rscale - 2.*R*b*Rscale*Rscale*cos((phiR - phib)*2.*M_PI),0.)), Qs)*2*M_PI*2*M_PI*kscale*Rscale*2*M_PI*Rscale*2*M_PI*pscale*(p*pscale+lambda); // bscale = Rscale //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
  return 0;
}


// Main program
int main(int argc, char *argv[]) {
  int rank;
  int size;

  int nev = 1;
  if (argc == 3) {
    nev = atoi(argv[2]);
  }

  // initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
  MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes

  int h5Flag = 0;
  pretty_ostream messager;

  display_logo();
  messager.flush("info");
  
  //printf("%.17f \n", TA(0.));

  // Cuba's parameters for integration
  int NDIM = 9;
  int NCOMP = 1;
  int NVEC = 1;
  //  double EPSREL = 5e-4;
  double EPSREL = 5e-3;
  double EPSABS = 1e-12;
  int VERBOSE = 0;
  int LAST = 4;
  int MINEVAL = 0;
  //int MAXEVAL = 100000000;
  int MAXEVAL = 100000000;
  int KEY = 0;
  
  //vegas
  int SEED = time(NULL);
  int NSTART = 10000000;
  int NINCREASE = 1000000;
  int NBATCH = 10000;
  int GRIDNO = 2;

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
 
  int comp, nregions, neval, fail;
  
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  params *userdata, data;

  data.Qs = 1.; // Saturation scale in GeV
  data.lambda = 0.1; // Infrared cutoff on p integral in GeV
  data.Y = 0.;
  userdata = &data; // Set the parameters to be passed to the integrand

  // // Run 9D Vegas integration
  // Vegas(NDIM, NCOMP, JPsiIntegrand1, userdata, NVEC,
  //       EPSREL, EPSABS, VERBOSE, SEED,
  //       MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
  //       GRIDNO, NULL, NULL,
  //       &neval, &fail, integral, error, prob);
  
  // // Print the result
  // double JPsi1result = (double)integral[0];
  // double JPsi1error = (double)error[0];
  // printf("%.8f +- %.8f\t\n", JPsi1result, JPsi1error);

  NDIM = 10;

  // Run 10D Vegas integration
  Vegas(NDIM, NCOMP, JPsiIntegrand2, userdata, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, NULL, NULL,
        &neval, &fail, integral, error, prob);
  
  // Print the result
  double JPsi2result = (double)integral[0];
  double JPsi2error = (double)error[0];
  printf("%.8f +- %.8f\t\n", JPsi2result, JPsi2error);




  NDIM = 8;

  // Run 8D Vegas integration
  Vegas(NDIM, NCOMP, FullIntegrand, userdata, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, NULL, NULL,
        &neval, &fail, integral, error, prob);
  
  // Print the result
  double gresult = (double)integral[0];
  double gerror = (double)error[0];
  printf("%.8f +- %.8f\t\n", gresult, gerror);
  

  // // Integrate 7D to get |p|-spectrum
  // int ppoints = 20; // Points in |p| to compute
  // double pstep = 0.1; // Step width in |p|
  
  // NDIM = 7;
  // int runs = 1;
  // for (int r=0; r<runs; r++){
  //   for (int i=1; i<=ppoints; i++){
  //     data.pe = i*pstep;
  //     userdata = &data;
  //     SEED = time(NULL)+r*10000;

  //     // Suave(NDIM, NCOMP, Integrand, userdata, NVEC,
  //     //       EPSREL, EPSABS, VERBOSE | LAST, SEED,
  //     //       MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
  //     //       NULL, NULL,
  //     //       &nregions, &neval, &fail, integral, error, prob);
      
  //     // Divonne(NDIM, NCOMP, Integrand, userdata, NVEC,
  //     //         EPSREL, EPSABS, VERBOSE, SEED,
  //     //         MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
  //     //         BORDER, MAXCHISQ, MINDEVIATION,
  //     //         NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
  //     //         NULL, NULL,
  //     //         &nregions, &neval, &fail, integral, error, prob);

  //     Vegas(NDIM, NCOMP, Integrand, userdata, NVEC,
  //           EPSREL, EPSABS, VERBOSE, SEED,
  //           MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
  //           GRIDNO, NULL, NULL,
  //           &neval, &fail, integral, error, prob);
      
  //     // Cuhre(NDIM, NCOMP, Integrand, userdata, NVEC,
  //     //       EPSREL, EPSABS, VERBOSE | LAST,
  //     //       MINEVAL, MAXEVAL, KEY,
  //     //       NULL, NULL,
  //     //       &nregions, &neval, &fail, integral, error, prob);
      
  //     gresult = 2.*constants::alphas/constants::CF/data.pe/data.pe*(double)integral[0];
  //     gerror = 2.*constants::alphas/constants::CF/data.pe/data.pe*(double)error[0];
  //     printf("%.3f \t \t%.8f \t%.8f\n", data.pe, gresult, gerror);
  //   }
  // }
    

  MPI_Finalize();

  return 1;
}

void display_logo() {
  cout << endl;
  cout << "- compute JPsi production with fluctuations -------------------------------------------------------------------" << endl;
  cout << endl;
}

