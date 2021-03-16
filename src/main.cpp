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
}

// Parameters that need to be passed to the integrand
struct params 
{
  double pe;
  double k;
  double Qs;
  double lambda;
};

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
  
  return result/67.09678472225216694;// normalization for above parameters RA=6.6fm and d=0.546fm - adjust if parameters change;
}

// Unintegrated gluon distribution for the proton
double phip(double k, double R, double Qs){
  return k*k*constants::CF/(pow(2.*M_PI,3))/constants::alphas
    *2.*constants::CF*exp(R*R/(2*constants::Bp))*exp(-constants::CF
                                                     *exp(R*R/(2*constants::Bp))*k*k/(constants::CA*Qs*Qs))/(constants::CA*Qs*Qs);
}

// Unintegrated gluon distribution for lead
double phit(double k, double R, double Qs){
  return k*k*constants::CF/(pow(2.*M_PI,3))/constants::alphas*2.*constants::CF*exp(-constants::CF*k*k/(constants::CA*Qs*Qs*TA(R)))/(constants::CA*Qs*Qs*TA(R));
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
#define f ff[0]

  double kscale = 100.;
  double Rscale = 12./constants::hbarc;
  double p = static_cast<params*>(userdata)->pe;
  double Qs = static_cast<params*>(userdata)->Qs;

  f = 2*M_PI*k*kscale*R*Rscale*b*Rscale*phip(k*kscale, R*Rscale, Qs)*phit(sqrt(p*p + k*k*kscale*kscale - 2.*p*k*kscale*cos((phi - phik)*2.*M_PI)), sqrt(max(R*Rscale*R*Rscale + b*b*Rscale*Rscale - 2.*R*b*Rscale*Rscale*cos((phiR - phib)*2.*M_PI),0.)), Qs)*2*M_PI*2*M_PI*kscale*Rscale*2*M_PI*Rscale*2*M_PI; // bscale = Rscale //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube)
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

  f = 2.*constants::alphas/constants::CF/(p*pscale+lambda)/(p*pscale+lambda)*2.*M_PI*k*kscale*R*Rscale*b*Rscale*phip(k*kscale, R*Rscale, Qs)*phit(sqrt((p*pscale+lambda)*(p*pscale+lambda) + k*k*kscale*kscale - 2.*(p*pscale+lambda)*k*kscale*cos((phi - phik)*2.*M_PI)), sqrt(max(R*Rscale*R*Rscale + b*b*Rscale*Rscale - 2.*R*b*Rscale*Rscale*cos((phiR - phib)*2.*M_PI),0.)), Qs)*2*M_PI*2*M_PI*kscale*Rscale*2*M_PI*Rscale*2*M_PI*pscale*(p*pscale+lambda); // bscale = Rscale //scaled phi (and dphi) to 2 pi phi etc. (as integral is always over unit cube) 
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
  int NDIM = 8;
  int NCOMP = 1;
  int NVEC = 1;
  double EPSREL = 5e-4;
  double EPSABS = 1e-12;
  int VERBOSE = 0;
  int LAST = 4;
  int MINEVAL = 0;
  int MAXEVAL = 100000000;
  //int MAXEVAL = 1000000;
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

  userdata = &data; // Set the parameters to be passed to the integrand

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
  
  // Integrate 7D to get |p|-spectrum
  int ppoints = 20; // Points in |p| to compute
  double pstep = 0.1; // Step width in |p|
  
  NDIM = 7;
  int runs = 1;
  for (int r=0; r<runs; r++){
    for (int i=1; i<=ppoints; i++){
      data.pe = i*pstep;
      userdata = &data;
      SEED = time(NULL)+r*10000;

      // Suave(NDIM, NCOMP, Integrand, userdata, NVEC,
      //       EPSREL, EPSABS, VERBOSE | LAST, SEED,
      //       MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
      //       NULL, NULL,
      //       &nregions, &neval, &fail, integral, error, prob);
      
      // Divonne(NDIM, NCOMP, Integrand, userdata, NVEC,
      //         EPSREL, EPSABS, VERBOSE, SEED,
      //         MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
      //         BORDER, MAXCHISQ, MINDEVIATION,
      //         NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
      //         NULL, NULL,
      //         &nregions, &neval, &fail, integral, error, prob);

      Vegas(NDIM, NCOMP, Integrand, userdata, NVEC,
            EPSREL, EPSABS, VERBOSE, SEED,
            MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
            GRIDNO, NULL, NULL,
            &neval, &fail, integral, error, prob);
      
      // Cuhre(NDIM, NCOMP, Integrand, userdata, NVEC,
      //       EPSREL, EPSABS, VERBOSE | LAST,
      //       MINEVAL, MAXEVAL, KEY,
      //       NULL, NULL,
      //       &nregions, &neval, &fail, integral, error, prob);
      
      gresult = 2.*constants::alphas/constants::CF/data.pe/data.pe*(double)integral[0];
      gerror = 2.*constants::alphas/constants::CF/data.pe/data.pe*(double)error[0];
      printf("%.3f \t \t%.8f \t%.8f\n", data.pe, gresult, gerror);
    }
    
  }
  MPI_Finalize();

  return 1;
}

void display_logo() {
  cout << endl;
  cout << "- compute JPsi production with fluctuations -------------------------------------------------------------------" << endl;
  cout << endl;
}

