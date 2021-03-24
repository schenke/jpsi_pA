#include "TAInt.h"

// Constants to be used
namespace constants {
  const int Nc = 3;
  const double hbarc = 0.1973269804;
}

// Nuclear density for lead
double TAInt::rhoA(double z, void * params) {
  // do everything in GeV^some power
  double RA = 6.62/constants::hbarc;
  double d = 0.546/constants::hbarc;
  double R = *(double *) params;
  double f = 1/(1 + exp((sqrt(R*R + z*z) - RA)/d))/1.;
  return f;
}


void TAInt::computeTAIntegral(){
  for(int i=0; i<200; i++){
    double R = double(i)/200.*20./constants::hbarc;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double result, error;
    gsl_function F;
    F.function = &TAInt::rhoA;
    F.params = &R;
    gsl_integration_qagi(&F, 1e-12, 1e-7, 1000, w, &result, &error);
    gsl_integration_workspace_free (w);
    xgrid[i] = R;
    TAgrid[i] = result/67.09678472225216694;// normalization for above parameters RA=6.62fm and d=0.546fm - adjust if parameters change;
    //    cout << xgrid[i] << " " << TAgrid[i] << endl;
  }
}
  
double TAInt::returnTA(double R){
  // gsl_interp_accel *acc
  //   = gsl_interp_accel_alloc ();
  // gsl_spline *spline
  //   = gsl_spline_alloc (gsl_interp_cspline, 200);

  // gsl_spline_init (spline, xgrid, TAgrid, 200);

  // double TA = gsl_spline_eval (spline, R, acc);
  // gsl_spline_free (spline);
  // gsl_interp_accel_free (acc);

  // simple linear interpolation
  int i = int(R*200./20.*constants::hbarc);
  //cout << "i=" << i<< " R=" << R << endl;
  double TA = (i+1-R*200./20.*constants::hbarc)*TAgrid[i] + (R*200./20.*constants::hbarc-i)*TAgrid[i+1];
  if(R>20./constants::hbarc){
    //cout << "range problem: " << TA << endl;
    TA=0.;}
  return TA;
}



