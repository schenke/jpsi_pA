#include "MV.h"

// Constants to be used
namespace constants {
  const double PI = 3.14159265358979323846;
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
  double lambda = 0.2; // IR regulator in GeV
  double A = ((double *)params)[0];
  double k = ((double *)params)[1];
  double f = 2.*constants::PI*pow(2.718281828+(1./lambda/lambda)/z/z,-A*z*z)*z*gsl_sf_bessel_J0(z*k);
  //double f = 2.*constants::PI*exp(-A*z*z)*z*gsl_sf_bessel_J0(z*k); //GBW for testing the numerics only
  return f;
}


// Unintegrated gluon distribution in MV
void MV::computePhip(){
  double result, error;
  gsl_function F;
  F.function = &MV::MVintegrandForList;
  double sum=0.;
  int max_steps = 1000;
  double a,b;

  for (int iA=0; iA<sizeA; iA++){
    double A=iA*deltaA;
    //double A=exp(double(iA/150.))-1.;
    cout << iA << endl;
    for (int ik=0; ik<sizek; ik++){
      double k=ik*deltak+0.0001;
      sum=0.;
      double params[] = { A, k };
      F.params = params;  
#pragma omp parallel for private(result,a,b) reduction(+:sum)
     for(int n=0; n < max_steps; n++)
        {
         gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
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
                               1e-12,// eps absolute
                               1e-6,// eps relative
                               1000, GSL_INTEG_GAUSS15,
                               w,
                               &result,
                               &error);
          sum += result;
          gsl_integration_workspace_free(w);
        }
      
      if (sum<0.){
        sum=0.;
      }
      Phip_array[iA][ik] = sum;
    }
  }
   
}


double MV::Phip(double k, double R, double Qs){

  double A = constants::CA/4./constants::CF*exp(-R*R/2./constants::Bp)*Qs*Qs;
  
  //int iA = int(log(A+1.)*150.); 
  int iA = int(A/deltaA); 
  int ik = int((k+0.0001)/deltak);
  //  cout << iA << " " << ik << endl;
  
  if (iA>=sizeA){
    cerr << "MV::Phip: A out of range." << endl;
    return 0.;
  }
  if (ik>=sizek){
    cerr << "MV::Phip: k out of range." << endl;
    return 0.;
  }

  //  cout << iA << endl;

  double Phip1 = (double(iA+1)-(A/deltaA))*Phip_array[iA][ik] + (A/deltaA-double(iA))*Phip_array[iA+1][ik];
  double Phip2 = (double(iA+1)-(A/deltaA))*Phip_array[iA][ik+1] + (A/deltaA-double(iA))*Phip_array[iA+1][ik+1];
  double result = (double(ik+1)-((k+0.0001)/deltak))*Phip1 + ((k+0.0001)/deltak-double(ik))*Phip2;
  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}

double MV::Phit(double k, double TA, double Qs){

  double A = constants::CA/4./constants::CF*TA*Qs*Qs;
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.0001)/deltak);
  //  cout << iA << " " << ik << endl;
  
  if (iA>=sizeA){
    cerr << "MV::Phit: A out of range." << endl;
    return 0.;
  }
  if (ik>=sizek){
    cerr << "MV::Phit: k out of range for k = " << k << "." << endl;
    return 0.;
  }

  //  cout << iA << endl;

  double Phip1 = (double(iA+1)-(A/deltaA))*Phip_array[iA][ik] + (A/deltaA-double(iA))*Phip_array[iA+1][ik];
  double Phip2 = (double(iA+1)-(A/deltaA))*Phip_array[iA][ik+1] + (A/deltaA-double(iA))*Phip_array[iA+1][ik+1];
  double result = (double(ik+1)-((k+0.0001)/deltak))*Phip1 + ((k+0.0001)/deltak-double(ik))*Phip2;
  //double Phip1 = (double(iA+1)-(log(A+1.)*150.))*Phip_array[iA][ik] + (log(A+1.)*150.-double(iA))*Phip_array[iA+1][ik];
  //double Phip2 = (double(iA+1)-(log(A+1.)*150.))*Phip_array[iA][ik+1] + (log(A+1.)*150.-double(iA))*Phip_array[iA+1][ik+1];
  //double result = (double(ik+1)-((k+0.0001)/deltak))*Phip1 + ((k+0.0001)/deltak-double(ik))*Phip2;
  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}


double MV::StF(double k, double TA, double Qs){

  double A = Qs*Qs/4.*TA;
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.0001)/deltak);
  //  cout << iA << " " << ik << endl;
  
  if (iA>=sizeA){
    cerr << "MV::Stf: A out of range." << endl;
    return 0.;
  }
  if (ik>=sizek){
    cerr << "MV::Stf: k out of range for k = " << k << "." << endl;
    return 0.;
  }

  //  cout << iA << endl;
  double Phip1 = (double(iA+1)-(A/deltaA))*Phip_array[iA][ik] + (A/deltaA-double(iA))*Phip_array[iA+1][ik];
  double Phip2 = (double(iA+1)-(A/deltaA))*Phip_array[iA][ik+1] + (A/deltaA-double(iA))*Phip_array[iA+1][ik+1];
  double result = (double(ik+1)-((k+0.0001)/deltak))*Phip1 + ((k+0.0001)/deltak-double(ik))*Phip2;

  // double Phip1 = (double(iA+1)-(log(A+1.)*150.))*Phip_array[iA][ik] + (log(A+1.)*150.-double(iA))*Phip_array[iA+1][ik];
  // double Phip2 = (double(iA+1)-(log(A+1.)*150.))*Phip_array[iA][ik+1] + (log(A+1.)*150.-double(iA))*Phip_array[iA+1][ik+1];
  // double result = (double(ik+1)-((k+0.0001)/deltak))*Phip1 + ((k+0.0001)/deltak-double(ik))*Phip2;
  //cout << "Phip=" << result << endl;
  return result;
}


int MV::writeTable(){
  
  stringstream name;
  name << "MVTable.dat";
  string MVTableName;
  MVTableName = name.str();
  
  ofstream Outfile1;
  Outfile1.open(MVTableName.c_str(), ios::out | ios::binary);

  sizeA = 800;
  sizek = 400;
  deltaA = 1./80.;
  deltak = 1./10.;
  
  // print header ------------- //
  Outfile1.write((char *)&sizeA, sizeof(int));
  Outfile1.write((char *)&sizek, sizeof(int));
  Outfile1.write((char *)&deltaA, sizeof(double));
  Outfile1.write((char *)&deltak, sizeof(double));
  
  double *val1 = new double[1];

  for (int iA=0; iA<sizeA; iA++){
    double A=iA*deltaA;
    for (int ik=0; ik<sizek; ik++){
      double k=ik*deltak+0.0001;
      val1[0] = Phip_array[iA][ik];
      Outfile1.write((char *)val1, sizeof(double));
    }
  }
  
  delete val1;
  
  if (Outfile1.good() == false) {
    std::cerr << "Error -- binary output of MV Table failed."
              << std::endl;
    return 0;
  }
  
  Outfile1.close();

  cout << "Wrote " << MVTableName
       << endl;
  return 1;
}


int MV::readTable(){
  stringstream name;
  name << "MVTable.dat";
  string MVTableName;
  MVTableName = name.str();

  std::ifstream InStream;
  InStream.precision(10);
  InStream.open(MVTableName, std::ios::in | std::ios::binary);
  cout << "Reading MV table from file " << MVTableName << endl;
  
  if(InStream.is_open())
    {
      // read parameters
      double temp;
      InStream.read(reinterpret_cast<char*>(&sizeA), sizeof(int));
      InStream.read(reinterpret_cast<char*>(&sizek), sizeof(int));
      InStream.read(reinterpret_cast<char*>(&deltaA), sizeof(double));
      InStream.read(reinterpret_cast<char*>(&deltak), sizeof(double));
      
      cout << "sizeA = " << sizeA << " sizek = " << sizek << " deltaA = " << deltaA << " deltak = " << deltak << endl;                               
                
      // read data
      double ValueBuffer;
      int INPUT_CTR=0;
                
      while( InStream.read(reinterpret_cast<char*>(&ValueBuffer), sizeof(double)))
      {
        
        int PositionIndx = (INPUT_CTR);
        
        int iA = PositionIndx / sizek;
        int ik = PositionIndx - sizek*iA;
                
        Phip_array[iA][ik] = ValueBuffer;
        
        INPUT_CTR++;
      }
    }
  
  else
    {
      std::cerr << "Error. Could not open file " << MVTableName << endl;
      return 0;
    }
  InStream.close();
  
  return 1;
 }
