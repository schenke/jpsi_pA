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
  double lambda = 0.2; // IR regulator in GeV
  double A = ((double *)params)[0];
  double k = ((double *)params)[1];
  double f = 2.*constants::PI*pow(2.718281828+(1./lambda/lambda)/z/z,-A*z*z)*z*gsl_sf_bessel_J0(z*k);
  //double f = 2.*constants::PI*pow(2.718281828+(1./lambda/z),-A*z*z)*z*gsl_sf_bessel_J0(z*k);
  //double f = 2.*constants::PI*exp(-A*z*z)*z*gsl_sf_bessel_J0(z*k); //GBW for testing the numerics only
  return f;
}


// fakeBK r integrand (parametrization that has faster x evolution at small r than at large r) y=-log(x)
double MV::BKintegrandForList(double z, void * params) {
  double lambda = 0.2; // IR regulator in GeV
  double A = ((double *)params)[0];
  double k = ((double *)params)[1];
  double y = ((double *)params)[2];  
  double f = 2.*constants::PI*pow(2.718281828+(1./lambda)/z,-A*z*z*pow(0.322*pow(0.01/(exp(-y)),(0.3/2.*0.5*(2.*exp(-0.16*z)/pow(z,0.05)))),2.))*z*gsl_sf_bessel_J0(z*k);
  //double f = 2.*constants::PI*pow(2.718281828+(1./lambda/z),-A*z*z)*z*gsl_sf_bessel_J0(z*k);
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
      double k=ik*deltak+0.01;
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


// Unintegrated gluon distribution in fake BK
void MV::computePhipBK(){
  double result, error;
  gsl_function F;
  F.function = &MV::BKintegrandForList;
  double sum=0.;
  int max_steps = 1000;
  double a,b;

  for (int iA=0; iA<sizeA; iA++){
    double A=iA*deltaA;
    //double A=exp(double(iA/150.))-1.;
    cout << iA << endl;
    for (int ik=0; ik<sizek; ik++){
      double k=ik*deltak+0.01;
      for (int iy=0; iy<sizey; iy++){
        double y = double(iy)*deltay;
        sum=0.;
        double params[] = { A, k, y };
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
        Phip_arrayBK[iA][ik][iy] = sum;
      }
    }
  }
}



double MV::PhipFluc(double k, double Tp, double Qs, double sizeFactor){

  double A = constants::CA/4./constants::CF*Tp*Qs*Qs;
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.01)/deltak);
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
  double result = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;
  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}


double MV::PhipBKFluc(double k, double Tp, double x){

  double A = constants::CA/4./constants::CF*Tp*2.19; //2.19 goes from b-independent proton to Gaussian proto (normalizes to TA=1 at b=0)
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.01)/deltak);
  double y = -log(x);
  int iy = int(y);
  
  
  if (iA>=sizeA){
    cerr << "MV::PhipBK: A out of range." << endl;
    return 0.;
  }
  if (ik>=sizek){
    cerr << "MV::PhipBK: k out of range." << endl;
    return 0.;
  }
  if (iy>=sizey){
    cerr << "MV::PhipBK: x out of range." << endl;
    return 0.;
  }

  double Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy];
  double Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy];
  double result1 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy+1];
  Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy+1];
  double result2 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  double result = (double(iy+1)-y/deltay)*result1 + (y/deltay-double(iy))*result2;

  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}

double MV::PhipBK(double k, double R, double sizeFactor, double x){
  
  //double bfactor = 2.19; //2.19 goes from b-independent proton to Gaussian proto (normalizes to TA=1 at b=0)
  double bfactor = 1.;
  double A = constants::CA/4./constants::CF*exp(-R*R/2./(constants::Bp*sizeFactor))*bfactor; 
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.01)/deltak);
  double y = -log(x);
  int iy = int(y);
  
  if (iA>=sizeA){
    cerr << "MV::PhipBK: A out of range." << endl;
    return 0.;
  }
  if (ik>=sizek){
    cerr << "MV::PhipBK: k out of range." << endl;
    return 0.;
  }
  if (iy>=sizey){
    cerr << "MV::PhipBK: x out of range." << endl;
    return 0.;
  }

  double Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy];
  double Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy];
  double result1 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy+1];
  Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy+1];
  double result2 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  double result = (double(iy+1)-y/deltay)*result1 + (y/deltay-double(iy))*result2;

  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}


double MV::Phip(double k, double R, double Qs, double sizeFactor){

  double A = constants::CA/4./constants::CF*exp(-R*R/2./(constants::Bp*sizeFactor))*Qs*Qs;
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.01)/deltak);
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
  double result = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;
  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}



double MV::Phit(double k, double TA, double Qs){

  double A = constants::CA/4./constants::CF*TA*Qs*Qs;
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.01)/deltak);
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
  double result = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;
  return k*k*constants::Nc/4./constants::alphas*result;
}


double MV::PhitBK(double k, double TA, double x){
  
 // double bfactor = 2.19; //2.19 goes from b-independent proton to Gaussian proton (normalizes to TA=1 at b=0)
  double bfactor = 1.;
  double A = constants::CA/4./constants::CF*TA*bfactor; 
  int iA = int(A/deltaA); 
  int ik = int((k+0.01)/deltak);
  double y = -log(x);
  int iy = int(y);
  
  if (iA>=sizeA){
    cerr << "MV::PhitBK: A out of range." << endl;
    return 0.;
  }
  if (ik>=sizek){
    cerr << "MV::PhitBK: k out of range." << endl;
    return 0.;
  }
  if (iy>=sizey){
    cerr << "MV::PhitBK: x out of range." << endl;
    return 0.;
  }

  double Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy];
  double Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy];
  double result1 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy+1];
  Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy+1];
  double result2 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  double result = (double(iy+1)-y/deltay)*result1 + (y/deltay-double(iy))*result2;

  //cout << "Phip=" << result << endl;
  return k*k*constants::Nc/4./constants::alphas*result;
}


double MV::StF(double k, double TA, double Qs){

  double A = Qs*Qs/4.*TA;
  
  int iA = int(A/deltaA);
  int ik = int((k+0.01)/deltak);
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
  double result = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;
  return result;
}




double MV::StFBK(double k, double TA, double x){

  double A = 1./4.*TA*2.19;//2.19 goes from b-independent proton to Gaussian proto (normalizes to TA=1 at b=0); 
  
  int iA = int(A/deltaA); 
  int ik = int((k+0.01)/deltak);
  double y = -log(x);
  int iy = int(y);
  
  if (iA>=sizeA){
    cerr << "MV::StFBK: A out of range." << endl;
    return 0.;
  }
  if (ik>=sizek){
    cerr << "MV::StFBK: k out of range." << endl;
    return 0.;
  }
  if (iy>=sizey){
    cerr << "MV::StFBK: x out of range." << endl;
    return 0.;
  }

  double Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy];
  double Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy];
  double result1 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  Phip1 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik][iy+1];
  Phip2 = (double(iA+1)-(A/deltaA))*Phip_arrayBK[iA][ik+1][iy+1] + (A/deltaA-double(iA))*Phip_arrayBK[iA+1][ik+1][iy+1];
  double result2 = (double(ik+1)-((k+0.01)/deltak))*Phip1 + ((k+0.01)/deltak-double(ik))*Phip2;

  double result = (double(iy+1)-y/deltay)*result1 + (y/deltay-double(iy))*result2;

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
  
  // print header ------------- //
  Outfile1.write((char *)&sizeA, sizeof(int));
  Outfile1.write((char *)&sizek, sizeof(int));
  Outfile1.write((char *)&deltaA, sizeof(double));
  Outfile1.write((char *)&deltak, sizeof(double));
  
  double *val1 = new double[1];

  for (int iA=0; iA<sizeA; iA++){
    double A=iA*deltaA;
    for (int ik=0; ik<sizek; ik++){
      double k=ik*deltak+0.01;
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

int MV::writeTableBK(){
  
  stringstream name;
  name << "BKTable.dat";
  string BKTableName;
  BKTableName = name.str();
  
  ofstream Outfile1;
  Outfile1.open(BKTableName.c_str(), ios::out | ios::binary);
  
  // print header ------------- //
  Outfile1.write((char *)&sizeA, sizeof(int));
  Outfile1.write((char *)&sizek, sizeof(int));
  Outfile1.write((char *)&sizey, sizeof(int));
  Outfile1.write((char *)&deltaA, sizeof(double));
  Outfile1.write((char *)&deltak, sizeof(double));
  Outfile1.write((char *)&deltay, sizeof(double));
  
  double *val1 = new double[1];

  for (int iA=0; iA<sizeA; iA++){
    //    double A=iA*deltaA;
    for (int ik=0; ik<sizek; ik++){
      //double k=ik*deltak+0.01;
      for (int iy=0; iy<sizey; iy++){
        //double y=iy*deltay;
        val1[0] = Phip_arrayBK[iA][ik][iy];
        Outfile1.write((char *)val1, sizeof(double));
      }
    }
  }
  delete val1;
  
  if (Outfile1.good() == false) {
    std::cerr << "Error -- binary output of MV Table failed."
              << std::endl;
    return 0;
  }
  
  Outfile1.close();
  
  cout << "Wrote " << BKTableName
       << endl;
  return 1;
}



// int MV::writeTableText(){
  
//   stringstream name;
//   name << "MVTableTest.dat";
//   string MVTableName;
//   MVTableName = name.str();
  
//   ofstream Outfile1;
//   Outfile1.open(MVTableName.c_str(), ios::out);
  
//   // print header ------------- //
//   Outfile1 << (sizeA) << " ";
//   Outfile1 << (sizek) << " ";
//   Outfile1 << deltaA << " ";
//   Outfile1 << deltak << endl;
  
//   for (int iA=0; iA<sizeA; iA++){
//     double A=iA*deltaA;
//     for (int ik=0; ik<sizek; ik++){
//       double k=ik*deltak+0.01;
//       Outfile1 << Phip_array[iA][ik] << endl;
//     }
//   }
  
//   Outfile1.close();

//   cout << "Wrote " << MVTableName
//        << endl;
//   return 1;
// }



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

int MV::readTableBK(){
  stringstream name;
  name << "BKTable.dat";
  string BKTableName;
  BKTableName = name.str();

  std::ifstream InStream;
  InStream.precision(10);
  InStream.open(BKTableName, std::ios::in | std::ios::binary);
  cout << "Reading BK table from file " << BKTableName << endl;
  
  if(InStream.is_open())
    {
      // read parameters
      double temp;
      InStream.read(reinterpret_cast<char*>(&sizeA), sizeof(int));
      InStream.read(reinterpret_cast<char*>(&sizek), sizeof(int));
      InStream.read(reinterpret_cast<char*>(&sizey), sizeof(int));
      InStream.read(reinterpret_cast<char*>(&deltaA), sizeof(double));
      InStream.read(reinterpret_cast<char*>(&deltak), sizeof(double));
      InStream.read(reinterpret_cast<char*>(&deltay), sizeof(double));
      
      cout << "sizeA = " << sizeA << " sizek = " << sizek 
           << " sizey = " << sizey << " deltaA = " << deltaA 
           << " deltak = " << deltak << " deltay = " << deltay << endl;                               
                
      // read data
      double ValueBuffer;
      int INPUT_CTR=0;
                
      while( InStream.read(reinterpret_cast<char*>(&ValueBuffer), sizeof(double)))
      {
        
        int PositionIndx = (INPUT_CTR);
        
        int iA = PositionIndx / (sizek*sizey);
        int ik = (PositionIndx - (sizek*sizey)*iA)/sizey;
        int iy = PositionIndx - (sizek*sizey)*iA - ik*sizey;
       
        Phip_arrayBK[iA][ik][iy] = ValueBuffer;
        
        INPUT_CTR++;
      }
    }
  
  else
    {
      std::cerr << "Error. Could not open file " << BKTableName << endl;
      return 0;
    }
  InStream.close();
  
  return 1;
 }
