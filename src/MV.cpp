#include "MV.h"

// Constants to be used
namespace constants {
  const int Nc = 3;
  const double hbarc = 0.1973269804;
  const double CA = double(Nc);
  const double CF = (double(Nc)*double(Nc) - 1.)/(2.*double(Nc));
  const double alphas = 0.3;
  const double Bp = 4.;
}


// MV r integrand
double MV::MVintegrandForList(double z, void * params) {
  double CA = constants::CA;
  double CF = constants::CF;
  double Bp = constants::Bp;
  double A = ((double *)params)[0];
  double k = ((double *)params)[1];
  double f = 2.*M_PI*pow((2.718281828+25./z/z),-A*z*z)*z*gsl_sf_bessel_J0(z*k);
  return f;
}

// Unintegrated gluon distribution for the proton in MV
void MV::computePhip(){
  cout << "Generating unintegrated gluon distribution for the proton" << endl;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result, error;
  gsl_function F;
  F.function = &MV::MVintegrandForList;
  double sum=0.;
  int max_steps = 1000;
  double a,b;

  for (int iA=0; iA<200; iA++){
    double A=iA/20.;
    cout << iA << endl;
    for (int ik=0; ik<200; ik++){
      double k=ik/10.+0.01;
      
      sum=0.;
      double params[] = { A, k };
      F.params = params;  
      for(int n=0; n < max_steps; n++)
        {
          if(n==0){
            a=0.;
            b=gsl_sf_bessel_zero_J0(1)/k;
          }
          else{
            a=gsl_sf_bessel_zero_J0(n)/k;
            b=gsl_sf_bessel_zero_J0(n+1)/k;
          }
          
          gsl_integration_qag (&F,  // function
                               a,   // from
                               b,   // to
                               0,// eps absolute
                               1e-6,// eps relative
                               1000, GSL_INTEG_GAUSS15,
                               w,
                               &result,
                               &error);
          sum += result;
        }
      
      if (sum<0.){
        sum=0.;
      }
      
      //cout << double(iA)/20. << " " << double(ik)/20. << " " << sum << endl;
      Phip_array[iA][ik] = sum;
    }
  }
   
  gsl_integration_workspace_free (w);
}

double MV::Phip(double k, double R, double Qs){

  double A = constants::CA/4./constants::CF*exp(-R*R/2./constants::Bp)*Qs*Qs;
  
  int iA = int(A*20.);
  int ik = int((k+0.01)*10.);
  //  cout << iA << " " << ik << endl;
  
  if (iA>=200){
    cerr << "MV::Phip: A out of range." << endl;
    return 0.;
  }
  if (ik>=200){
    cerr << "MV::Phip: k out of range." << endl;
    return 0.;
  }

  //  cout << iA << endl;

  double Phip1 = (double(iA+1)-(A*20.))*Phip_array[iA][ik] + (A*20-double(iA))*Phip_array[iA+1][ik];
  double Phip2 = (double(iA+1)-(A*20.))*Phip_array[iA][ik+1] + (A*20.-double(iA))*Phip_array[iA+1][ik+1];
  double result = (double(ik+1)-((k+0.01)*10.))*Phip1 + ((k+0.01)*10.-double(ik))*Phip2;
  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}






