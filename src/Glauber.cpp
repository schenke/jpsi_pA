#include "Glauber.h"

using namespace std;

Glauber::Glauber(Parameters *inParam)
{
  param = new Parameters();
  param = inParam; 
  gslRan = gsl_rng_alloc (gsl_rng_taus);
  numberOfQuarks = 3;
}

// destructor
Glauber::~Glauber()
{
  delete param;
  delete hadronBins;
}

double Glauber::FermiDistribution(Nucleus *nuc, double r)
{
 double f;

 f = r*r / (1. + exp( (r-nuc->R_WS) / nuc->a_WS) );
  
 return f;
}

double Glauber::ExponentialDistribution(double a, double r)
{
 double f;
 //a = \sqrt{12}/R_p = 3.87, with R_p = 0.895 from http://journals.aps.org/rmp/pdf/10.1103/RevModPhys.77.1 
 f = r*r * exp(-a*r);
  
 return f;
}


// this is only for testing to see if the probability distributions are sampled correctly
void Glauber::makeHistograms(Random *random)
{
  double rbins[100];
  double rqbins[100];
  double xbins[100];
  double ybins[100];
  double zbins[100];
  double thetabins[40];
  double Rmax=15.;
  double x,y,z,r,rq, theta;
  int samples = 100000;

  Quark quark; 
  
  for(int j=0; j<99; j++)
    {
      rbins[j]=0.;
      rqbins[j]=0.;
      xbins[j]=0.;
      ybins[j]=0.;
      zbins[j]=0.;
    }
  for(int j=0; j<40; j++)
    {
      thetabins[j]=0.;
    }
    
  for(int i=0; i<samples; i++)
    {
      quark = sampleQuark(random);
      Target.nucleonList.push_back(sampleRho(&Target, random));
      // cout << "x=" << Target.nucleonList.at(i).x << ", y=" << Target.nucleonList.at(i).y << ", z=" << Target.nucleonList.at(i).z << endl;
      x=Target.nucleonList.at(i).x;
      y=Target.nucleonList.at(i).y;
      z=Target.nucleonList.at(i).z;
      r = sqrt(x*x+y*y+z*z);
      rq = sqrt(quark.x*quark.x+quark.y*quark.y+quark.z*quark.z);
      theta = atan2(y,x);
      //cout << r<< endl;
      for(int j=0; j<99; j++)
	{
	  if(r<Rmax/100.*(j+1) && r>Rmax/100.*j)
	    rbins[j]+=1./samples/(Rmax/100.);
	  if(rq<Rmax/100.*(j+1) && rq>Rmax/100.*j)
	    rqbins[j]+=1./samples/(Rmax/100.);
	  if(x<Rmax/100.*(j+1)-Rmax/2. && x>Rmax/100.*j-Rmax/2.)
	    xbins[j]+=1./samples/(Rmax/100.);
	  if(y<Rmax/100.*(j+1)-Rmax/2. && y>Rmax/100.*j-Rmax/2.)
	    ybins[j]+=1./samples/(Rmax/100.);
	  if(z<Rmax/100.*(j+1)-Rmax/2. && z>Rmax/100.*j-Rmax/2.)
	    zbins[j]+=1./samples/(Rmax/100.);
	}
      for(int j=0; j<39; j++)
	{
	  if(theta<M_PI/40.*(j+1) && theta>M_PI/40.*j)
	    thetabins[j]+=1./samples/(M_PI/40.);
	}
    }

    
  ofstream hist;
  stringstream ss_hist_name;
  ss_hist_name << "rDist.dat";
  string hist_name = ss_hist_name.str();
  hist.open(hist_name.c_str()); 
  
  for(int j=0; j<99; j++)
    {
      hist << Rmax/100*(j+0.5) << " " << rbins[j] << endl;
    }
  
  hist.close();

  ofstream qhist;
  stringstream ss_qhist_name;
  ss_qhist_name << "quarkDist.dat";
  string qhist_name = ss_qhist_name.str();
  qhist.open(qhist_name.c_str()); 
  
  for(int j=0; j<99; j++)
    {
      qhist << Rmax/100*(j+0.5) << " " << rqbins[j] << endl;
    }
  
  qhist.close();

  ofstream xhist;
  stringstream ss_xhist_name;
  ss_xhist_name << "xDist.dat";
  string xhist_name = ss_xhist_name.str();
  xhist.open(xhist_name.c_str()); 
  
  for(int j=0; j<99; j++)
    {
      xhist << Rmax/100*(j+0.5)-Rmax/2. << " " << xbins[j]  << " " << ybins[j]  << " " << zbins[j] << endl;
    }
  
  xhist.close();

  ofstream thetahist;
  stringstream ss_thetahist_name;
  ss_thetahist_name << "thetaDist.dat";
  string thetahist_name = ss_thetahist_name.str();
  thetahist.open(thetahist_name.c_str()); 
  
  for(int j=0; j<39; j++)
    {
      thetahist << M_PI/40.*(j+0.5) << " " << thetabins[j] << endl;
    }
  
  thetahist.close();
}




void Glauber::init(Random *random)
{
  if(!setNucleusParameters(&Target, param->getTarget()))
    {
      cerr << "Could not read Nucleus information for target " << param->getTarget() << endl;
      exit(1);
    }
     
  if(!setNucleusParameters(&Projectile, param->getProjectile()))
    {
      cerr << "Could not read Nucleus information for projectile " << param->getProjectile() << endl;
      exit(1);
    }

  int rank = MPI::COMM_WORLD.Get_rank(); //number of current processor
} 

void Glauber::makeNuclei(Random *random)
{
  int Nx = param->getOutputNumberOfTransverseCells();
  //reset hadron distribution
  for(int j=0; j<Nx*Nx; j++)
    {

      hadronBins[j] = 0.;
    }

  Target.nucleonList.clear();
  Projectile.nucleonList.clear();
  
  if( Target.A == 1 )
    {
      Nucleon nucleon;
      nucleon.x = 0.;
      nucleon.y = 0.;
      nucleon.z = 0.;
      nucleon.collided = 0;
      Target.nucleonList.push_back(nucleon);
    }
  else
    {
      for(int i=0; i<Target.A; i++)
	{
	  Target.nucleonList.push_back(sampleRho(&Target, random));
	}
      
      //shift to center of mass
      double meanx=0, meany=0, meanz=0;
      
      for(int i=0; i<Target.A; i++)
	{
	  meanx += Target.nucleonList.at(i).x;
	  meany += Target.nucleonList.at(i).y;
	  meanz += Target.nucleonList.at(i).z;
	}
      
      meanx /= double(Target.A);
      meany /= double(Target.A);
      meanz /= double(Target.A);
      
      for(int i=0; i<Target.A; i++)
	{
	  Target.nucleonList.at(i).x -= meanx;
	  Target.nucleonList.at(i).y -= meany;
	  Target.nucleonList.at(i).z -= meanz;
	}
      
      // // test if shift worked:
      // meanx = 0;
      // meany = 0;
      // meanz = 0;

      // for(int i=0; i<Target.A; i++)
      // 	{
      // 	  meanx += Target.nucleonList.at(i).x;
      // 	  meany += Target.nucleonList.at(i).y;
      // 	  meanz += Target.nucleonList.at(i).z;
      // 	}

      // for(int i=0; i<Target.A; i++)
      // 	{
      // 	  cout << "mean x = " << meanx << endl;
      // 	  cout << "mean y = " << meany << endl;
      // 	  cout << "mean z = " << meanz << endl;
      // 	}
           

    }
  if( Projectile.A == 1 )
    {
      Nucleon nucleon;
      nucleon.x = 0.;
      nucleon.y = 0.;
      nucleon.z = 0.;
      nucleon.collided = 0;
      Projectile.nucleonList.push_back(nucleon);
    }
  else
    {
      for(int i=0; i<Projectile.A; i++)
	{
	  Projectile.nucleonList.push_back(sampleRho(&Projectile, random));
	}
      
      //shift to center of mass
      double meanx=0, meany=0, meanz=0;
      
      for(int i=0; i<Projectile.A; i++)
	{
	  meanx += Projectile.nucleonList.at(i).x;
	  meany += Projectile.nucleonList.at(i).y;
	  meanz += Projectile.nucleonList.at(i).z;
	}
      
      meanx /= double(Projectile.A);
      meany /= double(Projectile.A);
      meanz /= double(Projectile.A);
      
      for(int i=0; i<Projectile.A; i++)
	{
	  Projectile.nucleonList.at(i).x -= meanx;
	  Projectile.nucleonList.at(i).y -= meany;
	  Projectile.nucleonList.at(i).z -= meanz;
	}
      
     
      // // test if shift worked:
      // meanx = 0;
      // meany = 0;
      // meanz = 0;

      // for(int i=0; i<Projectile.A; i++)
      // 	{
      // 	  meanx += Projectile.nucleonList.at(i).x;
      // 	  meany += Projectile.nucleonList.at(i).y;
      // 	  meanz += Projectile.nucleonList.at(i).z;
      // 	}

      // for(int i=0; i<Projectile.A; i++)
      // 	{
      // 	  cout << "mean x = " << meanx << endl;
      // 	  cout << "mean y = " << meany << endl;
      // 	  cout << "mean z = " << meanz << endl;
      // 	}
      
    }



  // if using quarks determine the quark positions within each nucleon
  if(param->getUseQuarks())
    {
      double xQuark[numberOfQuarks];
      if( Target.A == 1 )
	{
	  for(int j=0; j<numberOfQuarks; j++)
	    {
	      Target.nucleonList.at(0).quarkList.push_back(sampleQuark(random));
	    }
	}
      else
	{
	  for(int i=0; i<Target.A; i++)
	    {
	      for(int j=0; j<numberOfQuarks; j++)
		{
		  Target.nucleonList.at(i).quarkList.push_back(sampleQuark(random));
		}
	    }
	}

      if( Projectile.A == 1 )
	{
	  for(int j=0; j<numberOfQuarks; j++)
	    {
	      Projectile.nucleonList.at(0).quarkList.push_back(sampleQuark(random));
	    }
	}
      else
	{
	  for(int i=0; i<Projectile.A; i++)
	    {
	      for(int j=0; j<numberOfQuarks; j++)
		{
		  Projectile.nucleonList.at(i).quarkList.push_back(sampleQuark(random));
		}
	    }
	}
    } 
   
  for (int i = 0; i<Target.A; i++) // shift the target's position by -b/2
    {
      Target.nucleonList.at(i).x=Target.nucleonList.at(i).x-param->getb()/2.;
    }   

  for (int i = 0; i<Projectile.A; i++) // shift the projectile's position by +b/2 
    {
      Projectile.nucleonList.at(i).x=Projectile.nucleonList.at(i).x+param->getb()/2.;
    }   
}


int Glauber::setNucleusParameters(Nucleus *nuc, string name) {
  string densityFunction;
  int success=1;
  if (name.compare("Au") == 0) {
    nuc->A = 197;
    nuc->Z = 79;
    densityFunction = "3Fermi";
    nuc->R_WS = 6.37;
    nuc->w_WS = 0;
    nuc->a_WS = 0.535;
  } else if (name.compare("Pb") == 0) {
    nuc->A = 208.;
    nuc->Z = 82.;
    densityFunction = "3Fermi";
    nuc->R_WS = 6.62;
    nuc->w_WS = 0.;
    nuc->a_WS = 0.546;
  } else if (name.compare("p") == 0) {
    nuc->A = 1.;
    nuc->Z = 1.;
    densityFunction = "3Fermi";
    nuc->R_WS = 1.;
    nuc->w_WS = 0.;
    nuc->a_WS = 1.;
  } else if (name.compare("He3") == 0) {
    nuc->A = 3;
    nuc->Z = 2;
    densityFunction = "readFromFile";
    nuc->R_WS = 0;
    nuc->w_WS = 0;
    nuc->a_WS = 0;
  } else if (name.compare("d") == 0) {
    nuc->A = 2;
    nuc->Z = 1;
    densityFunction = "Hulthen";
    nuc->R_WS = 1.0;
    nuc->w_WS = 1.18;
    nuc->a_WS = 0.228;
  } else if (name.compare("C") == 0) {
    nuc->A = 12;
    nuc->Z = 6;
    densityFunction = "2HO";
    nuc->R_WS = 2.44;
    nuc->w_WS = 1.403;
    nuc->a_WS = 1.635;
  } else if (name.compare("O") == 0) {
    nuc->A = 16;
    nuc->Z = 8;
    densityFunction = "3Fermi";
    nuc->R_WS = 2.608;
    nuc->w_WS = -0.051;
    nuc->a_WS = 0.513;
  } else if (name.compare("S") == 0) {
    nuc->A = 32;
    nuc->Z = 16;
    densityFunction = "3Gauss";
    nuc->R_WS = 2.54;
    nuc->w_WS = 0.16;
    nuc->a_WS = 2.191;
  } else if (name.compare("W") == 0) {
    nuc->A = 184;
    nuc->Z = 74;
    densityFunction = "3Fermi";
    nuc->R_WS = 6.51;
    nuc->w_WS = 0;
    nuc->a_WS = 0.535;
  } else if (name.compare("Al") == 0) {
    nuc->A = 27;
    nuc->Z = 13;
    densityFunction = "3Fermi";
    nuc->R_WS = 3.07;
    nuc->w_WS = 0;
    nuc->a_WS = 0.519;
  } else if (name.compare("Ca") == 0) {
    nuc->A = 40;
    nuc->Z = 20;
    densityFunction = "3Fermi";
    nuc->R_WS = 3.766;
    nuc->w_WS = -0.161;
    nuc->a_WS = 0.586;
  } else if (name.compare("Cu") == 0) {
    nuc->A = 63;
    nuc->Z = 29;
    densityFunction = "3Fermi";
    nuc->R_WS = 4.163;
    nuc->w_WS = 0;
    nuc->a_WS = 0.606;
  } else if (name.compare("Fe") == 0) {
    nuc->A = 56;
    nuc->Z = 26;
    densityFunction = "3Fermi";
    nuc->R_WS = 4.106;
    nuc->w_WS = 0;
    nuc->a_WS = 0.519;
  } else if (name.compare("Pt") == 0) {
    nuc->A = 195;
    nuc->Z = 78;
    densityFunction = "3Fermi";
    nuc->R_WS = 6.78;
    nuc->w_WS = 0;
    nuc->a_WS = 0.54;
  } else if (name.compare("U") == 0) {
    nuc->A = 238;
    nuc->Z = 92;
    densityFunction = "3Fermi";
    nuc->R_WS = 6.81;
    nuc->w_WS = 0;
    nuc->a_WS = 0.55;
  } else if (name.compare("Ru") == 0) {
    nuc->A = 96;
    nuc->Z = 44;
    densityFunction = "3Fermi";
    nuc->R_WS = 5.085;
    nuc->w_WS = 0;
    nuc->a_WS = 0.46;
  } else if (name.compare("Zr") == 0) {
    nuc->A = 96;
    nuc->Z = 40;
    densityFunction = "3Fermi";
    nuc->R_WS = 5.02;
    nuc->w_WS = 0;
    nuc->a_WS = 0.46;
  } else if (name.compare("Xe") == 0) {
    nuc->A = 129;
    nuc->Z = 54;
    densityFunction = "3Fermi";
    nuc->R_WS = 5.42; // 5.36       // new values from arXiv:1508.06294
    nuc->w_WS = 0;
    nuc->a_WS = 0.57;    // 0.590;     // new values from arXiv:1508.06294
  }
  else
    success = 0;


 if(densityFunction.compare("2HO")==0) 
    {
      nuc->DensityFunc = 1; //NuInt2HO;
    }
  else if(densityFunction.compare("3Gauss")==0) 
    {
      nuc->DensityFunc = 2; //NuInt3Gauss;
    }
  else if(densityFunction.compare("3Fermi")==0) 
    {
      nuc->DensityFunc = 3; //NuInt3Fermi;
    }
  else if(densityFunction.compare("Hulthen")==0) 
    {
      nuc->DensityFunc = 8; //NuIntHulthen;  
    }
  else if(densityFunction.compare("readFromFile")==0) 
    {
      nuc->DensityFunc = 1;  
    }

  return success;
}

// int Glauber::setNucleusParameters(Nucleus *nucleus, string name)
// {
//   int success = 0;
//   stringstream convert;
//   ifstream inputFile("nuclei");
//   string line, variable, value, densityFunction;
//   while ( getline (inputFile,line) )
//    {
//      if ( param->strWord(1, line) == "Name" )
//        if ( param->strWord(2, line) == name )
//        { 
// 	 success = 1;
//  	 nucleus->name = name;
//  	 getline (inputFile,line);
// 	 if (param->strWord(1, line) == "A") 
// 	   {
// 	     convert << param->strWord(2, line);
// 	     convert >> nucleus->A;
// 	   }
// 	 else
// 	   {
// 	     cerr << "Error reading A. in Glauber::setNucleusParameters. please check structure of nuclei table. Exiting." << endl;
// 	     exit(1);
// 	   }
// 	 getline (inputFile,line);
// 	 convert.str(string());
// 	 convert.clear();
// 	 if (param->strWord(1, line) == "Z") 
// 	   {
// 	     convert << param->strWord(2, line);
// 	     convert >> nucleus->Z;
// 	   }
// 	 else
// 	   {
// 	     cerr << "Error reading Z. in Glauber::setNucleusParameters. please check structure of nuclei table. Exiting." << endl;
// 	     exit(1);
// 	   }
// 	 getline (inputFile,line);
// 	 if (param->strWord(1, line) == "density_func") 
// 	   densityFunction = param->strWord(2, line);
// 	 else
// 	   {
// 	     cerr << "Error reading density_func. in Glauber::setNucleusParameters. please check structure of nuclei table. Exiting." << endl;
// 	     exit(1);
// 	   }
// 	 getline (inputFile,line);
// 	 convert.str(string());
// 	 convert.clear();
// 	 if (param->strWord(1, line) == "R_WS") 
// 	   {
// 	     convert << param->strWord(2, line);
// 	     convert >> nucleus->R_WS;
// 	   }
// 	 else
// 	   {
// 	     cerr << "Error reading R_WS. in Glauber::setNucleusParameters. please check structure of nuclei table. Exiting." << endl;
// 	     exit(1);
// 	   }
// 	 getline (inputFile,line);
// 	 convert.str(string());
// 	 convert.clear();
// 	 if (param->strWord(1, line) == "w_WS") 
// 	   {
// 	     convert << param->strWord(2, line);
// 	     convert >> nucleus->w_WS;
// 	   }
// 	 else
// 	   {
// 	     cerr << "Error reading w_WS. in Glauber::setNucleusParameters. please check structure of nuclei table. Exiting." << endl;
// 	     exit(1);
// 	   }
// 	 getline (inputFile,line);
// 	 convert.str(string());
// 	 convert.clear();
// 	 if (param->strWord(1, line) == "a_WS") 
// 	   {
// 	     convert << param->strWord(2, line);
// 	     convert >> nucleus->a_WS;
// 	   }
// 	 else
// 	   {
// 	     cerr << "Error reading a_WS. in Glauber::setNucleusParameters. please check structure of nuclei table. Exiting." << endl;
// 	     exit(1);
// 	   }
//        }
//    }
//   inputFile.close(); 

//   if(densityFunction.compare("2HO")==0) 
//     {
//       nucleus->DensityFunc = 1; //NuInt2HO;
//     }
//   else if(densityFunction.compare("3Gauss")==0) 
//     {
//       nucleus->DensityFunc = 2; //NuInt3Gauss;
//     }
//   else if(densityFunction.compare("3Fermi")==0) 
//     {
//       nucleus->DensityFunc = 3; //NuInt3Fermi;
//     }
//   else if(densityFunction.compare("Hulthen")==0) 
//     {
//       nucleus->DensityFunc = 8; //NuIntHulthen;  
//     }
//   else if(densityFunction.compare("readFromFile")==0) 
//     {
//       nucleus->DensityFunc = 1;  
//     }
  
//   if(success==0)
//     return 0;
//   else
//     return 1;
// }


Nucleon Glauber::sampleRho(Nucleus *nuc, Random *random)
{
  Nucleon nucleon;
  double r, x, y, z, tmp;
  double phi, theta;
  double a_WS = nuc->a_WS;
  double R_WS = nuc->R_WS;
   
  phi = 2.*M_PI*random->genrand64_real1();
  theta = acos(1. - 2.*random->genrand64_real1());
 
  do
    {
      // sample the radius from the envelope distribution
      r = 2.*a_WS*R_WS*sqrt(-log(-((-1 + random->genrand64_real1() * ((1.-exp(-(1000000/(4.*a_WS*a_WS*R_WS*R_WS)))))))));
     
      // sample uniform random number to decide whether to accept or reject the sampled one
      tmp = random->genrand64_real1();

      // warn if the envelope happens to go below the actual distriution (should never happen)
      if ( FermiDistribution(nuc,r) > (6.*r)/a_WS*exp(-r*r/(4.*R_WS*R_WS*a_WS*a_WS)) )
	cerr << "WARNING: rho>envelope: " << "rho=" << FermiDistribution(nuc,r) 
	     << ", envelope=" <<  (6.*r)/a_WS*exp(-r*r/(4.*R_WS*R_WS*a_WS*a_WS)) << endl;
      //repeat until tmp is smaller than the ratio p(y)/f(y)
    } while( tmp > FermiDistribution(nuc,r) / ((6.*r)/a_WS * exp(-r*r/(4.*R_WS*R_WS*a_WS*a_WS))) ); 

  // determine x,y,z coordinates of the nucleon
  x=r*sin(theta)*cos(phi);
  y=r*sin(theta)*sin(phi);
  z=r*cos(theta);
  
  // set values of nucleon
  nucleon.x=x;
  nucleon.y=y;
  nucleon.z=z;
  nucleon.collided=0;
  
  return nucleon; 
}

Quark Glauber::sampleQuark(Random *random)
{
  Quark quark;
  double r, x, y, z, tmp;
  double phi, theta;
  double a = 3.87;
  double xQuark[numberOfQuarks];

  phi = 2.*M_PI*random->genrand64_real1();
  theta = acos(1. - 2.*random->genrand64_real1());
 
  do
    {
      // sample the radius from the envelope distribution
      r = sqrt(3)*sqrt(-log(1. + (-1. + 1./exp(1000000/3)) * random->genrand64_real1()));
      
      // sample uniform random number to decide whether to accept or reject the sampled one
      tmp = random->genrand64_real1();
      
      // warn if the envelope happens to go below the actual distriution (should never happen)
      if ( ExponentialDistribution(a,r) > (r/2.) / a * exp(-(r*r)/3.) )
	cerr << "WARNING: rho>envelope: " << "rho=" << ExponentialDistribution(a,r)
	     << ", envelope=" <<  (r/2.) / a * exp(-(r*r)/3.) << endl;
      //repeat until tmp is smaller than the ratio p(y)/f(y)
    } while( tmp > ExponentialDistribution(a,r) / (  (r/2.) / a * exp(-(r*r)/3.) )); 

  // determine x,y,z coordinates of the quark (relative to the nucleon)
  x=r*sin(theta)*cos(phi);
  y=r*sin(theta)*sin(phi);
  z=r*cos(theta);
  
  // set values of nucleon
  quark.x=x;
  quark.y=y;
  quark.z=z;
  quark.collided=0;

  return quark; 
}


// make roots a parameter and compute rapidity. also make rapidity dependent on the x of the quarks

int Glauber::collide(Random *random)
{
  double dx, dy, dij; 
  int Ncoll=0;
   
  // ------------------------------------- use nucleons -------------------------------------

  if(param->getUseQuarks() == 0)
    {
      //      double d2 = param->getLexusLambda()*param->getSigmaNNinel()/(M_PI*10.); // in fm^2 
      double d2 = param->getSigmaNNinel()/(M_PI*10.); // in fm^2 
      if (param->getGaussianWounding() == 0)
	{
	  for (int i = 0; i<Projectile.A; i++) 
	    {
	      for (int j = 0 ; j<Target.A ;j++) 
		{
		  dx = Target.nucleonList.at(j).x-Projectile.nucleonList.at(i).x;
		  dy = Target.nucleonList.at(j).y-Projectile.nucleonList.at(i).y;
		  dij = dx*dx+dy*dy;
		  if (dij < d2) 
		    {
		      Ncoll++;
		      Target.nucleonList.at(j).collided=1;
		      //Target.nucleonList.at(j).rapidity= -param->getBeamRapidity();
		      //      Target.nucleonList.at(j).collidedWith.push_back(i);           // commented 2/3/2015
		      Projectile.nucleonList.at(i).collided=1;
		      //     Projectile.nucleonList.at(i).rapidity=param->getBeamRapidity();
		      Projectile.nucleonList.at(i).collidedWith.push_back(j);
		    }
		}
	    }
	}
      else
	{
	  double p;
	  double G=0.92;
	  double ran;
	  
	  for (int i = 0; i<Projectile.A; i++) 
	    {
	      for (int j = 0 ; j<Target.A ;j++) 
		{
		  dx = Target.nucleonList.at(j).x-Projectile.nucleonList.at(i).x;
		  dy = Target.nucleonList.at(j).y-Projectile.nucleonList.at(i).y;
		  dij = dx*dx+dy*dy;
		  
		  p = G * exp(-G*dij/d2); // Gaussian profile 
		  
		  ran = random->genrand64_real1();
		  
		  if (ran < p) 
		    {
		      Ncoll++;
		      Target.nucleonList.at(j).collided=1;
		      //Target.nucleonList.at(j).rapidity= - param->getBeamRapidity();
		      //Target.nucleonList.at(j).collidedWith.push_back(i); // commented 2/3/2015
		      Projectile.nucleonList.at(i).collided=1;
		      //Projectile.nucleonList.at(i).rapidity=param->getBeamRapidity();
		      Projectile.nucleonList.at(i).collidedWith.push_back(j);
		    }
		}
	    }
	}
      //      cout << "Ncoll = " << Ncoll << endl;
    }

  // ------------------------------------- use quarks -------------------------------------
  else 
    {
      double y;
      double xQuarkT[numberOfQuarks];
      double xQuarkP[numberOfQuarks];
      //      double d2 = param->getLexusLambda()*param->getSigmaQQinel()/(M_PI*10.); // in fm^2
      double d2 = param->getSigmaQQinel()/(M_PI*10.); // in fm^2
      
      if (param->getGaussianWounding() == 0)
	{
	  for (int i = 0; i<Projectile.A; i++) 
	    {
	      for (int j = 0 ; j<Target.A ;j++) 
		{
		  for (int iq = 0; iq<numberOfQuarks; iq++) 
		    {
		      for (int jq = 0 ; jq<numberOfQuarks ;jq++) 
			{
			  dx = (Target.nucleonList.at(j).x + Target.nucleonList.at(j).quarkList.at(jq).x )
			    - (Projectile.nucleonList.at(i).x + Projectile.nucleonList.at(i).quarkList.at(iq).x );
			  dy = (Target.nucleonList.at(j).y + Target.nucleonList.at(j).quarkList.at(jq).y )
			    - (Projectile.nucleonList.at(i).y + Projectile.nucleonList.at(i).quarkList.at(iq).y );
			  dij = dx*dx+dy*dy;
			  if (dij < d2) 
			    {
			      Ncoll++;
			      //note: in practice we only need to make one list - 2nd list is commented now:
			      // e.g. only keep projectile's information on what target nucleons it collided with and remove the other to save time
			      Target.nucleonList.at(j).quarkList.at(jq).collided = 1;
			      Target.nucleonList.at(j).collided = 1;
			      //  Target.nucleonList.at(j).rapidity = - param->getBeamRapidity();
			      // Target.nucleonList.at(j).collidedWith.push_back(i);  // commented 2/3/2015
			      // Target.nucleonList.at(j).quarkList.at(jq).collidedWith.push_back(i);     // commented 2/3/2015   // what nucleon did the quark collide with
			      //Target.nucleonList.at(j).quarkList.at(jq).collidedWithQuark.push_back(iq);  // commented 2/3/2015 // what quark in that nucleon did the quark coll.
			   
			      Projectile.nucleonList.at(i).quarkList.at(iq).collided = 1;
			      Projectile.nucleonList.at(i).collided = 1;
			      //Projectile.nucleonList.at(i).rapidity = param->getBeamRapidity();
			      Projectile.nucleonList.at(i).quarkList.at(iq).collidedWith.push_back(j);       // what nucleon did the quark collide with
			      Projectile.nucleonList.at(i).quarkList.at(iq).collidedWithQuark.push_back(jq); // what quark in that nucleon did the quark coll.	      
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  double p;
	  double G=0.92;
	  double ran;

	  for (int i = 0; i<Projectile.A; i++) 
	    {
	      for (int j = 0 ; j<Target.A ;j++) 
		{
		  for (int iq = 0; iq<numberOfQuarks; iq++) 
		    {
		      for (int jq = 0 ; jq<numberOfQuarks ;jq++) 
			{
			  dx = (Target.nucleonList.at(j).x + Target.nucleonList.at(j).quarkList.at(jq).x )
			    - (Projectile.nucleonList.at(i).x + Projectile.nucleonList.at(i).quarkList.at(iq).x );
			  dy = (Target.nucleonList.at(j).y + Target.nucleonList.at(j).quarkList.at(jq).y )
			    - (Projectile.nucleonList.at(i).y + Projectile.nucleonList.at(i).quarkList.at(iq).y );
			  dij = dx*dx+dy*dy;
			  
			  p = G * exp(-G*dij/d2); // Gaussian profile 
			  
			  ran = random->genrand64_real1();
		  
			  if (ran < p) 
			    {
			      Ncoll++;
			      //note: in practice we only need to make one list - 2nd one is commented now:
			      // e.g. only keep projectile's information on what target nucleons it collided with and remove the other to save time
			      Target.nucleonList.at(j).quarkList.at(jq).collided = 1;
			      Target.nucleonList.at(j).collided = 1;
			      //  Target.nucleonList.at(j).rapidity = - param->getBeamRapidity();
			      //Target.nucleonList.at(j).collidedWith.push_back(i); // commented 2/3/2015
			      //Target.nucleonList.at(j).quarkList.at(jq).collidedWith.push_back(i);  // commented 2/3/2015      // what nucleon did the quark collide with
			      //Target.nucleonList.at(j).quarkList.at(jq).collidedWithQuark.push_back(iq);  // commented 2/3/2015// what quark in that nucleon did the quark coll.
			   
			      Projectile.nucleonList.at(i).quarkList.at(iq).collided = 1;
			      Projectile.nucleonList.at(i).collided = 1;
			      //Projectile.nucleonList.at(i).rapidity = param->getBeamRapidity();
			      Projectile.nucleonList.at(i).quarkList.at(iq).collidedWith.push_back(j);       // what nucleon did the quark collide with
			      Projectile.nucleonList.at(i).quarkList.at(iq).collidedWithQuark.push_back(jq); // what quark in that nucleon did the quark coll.
			    }
			}
		    }
		}
	    }
	}
      //      cout << "Ncoll = " << Ncoll << endl;
    }// end use quarks
  //cout << "Ncoll = " << Ncoll << endl;
  return Ncoll;
}

int Glauber::getNpart()
{
  int Npart=0;
  //calculate N_part
  
  for (int i = 0; i<Projectile.A; i++) 
    {
      if(Projectile.nucleonList.at(i).collided)
	Npart++;
    }
  
  for (int j = 0 ; j<Target.A ;j++) 
    {
      if(Target.nucleonList.at(j).collided)
	Npart++;
    }
 
  //cout << "Npart = " << Npart << endl;
  return Npart;
   
}


int Glauber::getAverageNumberOfCollisions()
{
  int collisions=0;

  if(param->getUseQuarks() == 0)
    // ------------------------------------- use nucleons -------------------------------------
    {
      
      for (int i = 0; i<Projectile.A; i++) 
	{
	  collisions += Projectile.nucleonList.at(i).collidedWith.size();
	}
      
      // for (int j = 0 ; j<Target.A ;j++) 
      // 	{
      // 	  collisions += Target.nucleonList.at(j).collidedWith.size();
      // 	}
      
      //cout << "collisions = " << collisions << endl;
      // cout << "average number of collisions per nucleon: " << collisions << endl;
      //      return collisions/(Projectile.A+Target.A);
      return collisions/(Projectile.A);
    }
  else
    // -------------------------------------- use quarks --------------------------------------
    {
      
     for (int i = 0; i<Projectile.A; i++) 
	{
	  for (int iq = 0; iq<numberOfQuarks; iq++) 
	    {
	      collisions += Projectile.nucleonList.at(i).quarkList.at(iq).collidedWith.size();
	    }
	}
    
     // for (int j = 0 ; j<Target.A ;j++) 
     //   {
     // 	  for (int jq = 0; jq<3; jq++) 
     // 	    {
     // 	      collisions += Target.nucleonList.at(j).quarkList.at(jq).collidedWith.size();
     // 	    }
     //   }
     //cout << "average number of collisions per quark: " << collisions << endl;
     // return collisions/(Projectile.A+Target.A)/3.;
     return collisions/(Projectile.A)/double(numberOfQuarks);
    }
  
}



double Glauber::unitStep(double a)
{
  if (a>=0)
    return 1.;
  else
    return 0.;
}



void Glauber::outputNucleonPositions()
{
  ofstream positionsTW;
  stringstream ss_positionsTW_name;
  ss_positionsTW_name << "PositionsTargetWounded.dat";
  string positionsTW_name = ss_positionsTW_name.str();
  positionsTW.open(positionsTW_name.c_str()); 

  ofstream positionsTS;
  stringstream ss_positionsTS_name;
  ss_positionsTS_name << "PositionsTargetSpectator.dat";
  string positionsTS_name = ss_positionsTS_name.str();
  positionsTS.open(positionsTS_name.c_str()); 
  
  ofstream positionsPW;
  stringstream ss_positionsPW_name;
  ss_positionsPW_name << "PositionsProjectileWounded.dat";
  string positionsPW_name = ss_positionsPW_name.str();
  positionsPW.open(positionsPW_name.c_str()); 

  ofstream positionsPS;
  stringstream ss_positionsPS_name;
  ss_positionsPS_name << "PositionsProjectileSpectator.dat";
  string positionsPS_name = ss_positionsPS_name.str();
  positionsPS.open(positionsPS_name.c_str()); 
  
  for(int i=0; i<Target.A; i++)
    {
      if(Target.nucleonList.at(i).collided == 1)
	positionsTW << Target.nucleonList.at(i).x << " " << Target.nucleonList.at(i).y << " " << Target.nucleonList.at(i).z << endl;
      if(Target.nucleonList.at(i).collided == 0)
	positionsTS << Target.nucleonList.at(i).x << " " << Target.nucleonList.at(i).y << " " << Target.nucleonList.at(i).z << endl;
    }

 for(int i=0; i<Projectile.A; i++)
    {
      if(Projectile.nucleonList.at(i).collided == 1)
	positionsPW << Projectile.nucleonList.at(i).x << " " << Projectile.nucleonList.at(i).y << " " << Projectile.nucleonList.at(i).z << endl;
      if(Projectile.nucleonList.at(i).collided == 0)
	positionsPS << Projectile.nucleonList.at(i).x << " " << Projectile.nucleonList.at(i).y << " " << Projectile.nucleonList.at(i).z << endl;
    }

  positionsTW.close();
  positionsTS.close();
  positionsPW.close();
  positionsPS.close();

}

void Glauber::outputQuarkPositions()
{
  ofstream positionsQTW;
  stringstream ss_positionsQTW_name;
  ss_positionsQTW_name << "PositionsTargetWoundedQuarks.dat";
  string positionsQTW_name = ss_positionsQTW_name.str();
  positionsQTW.open(positionsQTW_name.c_str()); 

  ofstream positionsQTS;
  stringstream ss_positionsQTS_name;
  ss_positionsQTS_name << "PositionsTargetSpectatorQuarks.dat";
  string positionsQTS_name = ss_positionsQTS_name.str();
  positionsQTS.open(positionsQTS_name.c_str()); 
  
  ofstream positionsQPW;
  stringstream ss_positionsQPW_name;
  ss_positionsQPW_name << "PositionsProjectileWoundedQuarks.dat";
  string positionsQPW_name = ss_positionsQPW_name.str();
  positionsQPW.open(positionsQPW_name.c_str()); 

  ofstream positionsQPS;
  stringstream ss_positionsQPS_name;
  ss_positionsQPS_name << "PositionsProjectileSpectatorQuarks.dat";
  string positionsQPS_name = ss_positionsQPS_name.str();
  positionsQPS.open(positionsQPS_name.c_str()); 
  
  for(int i=0; i<Target.A; i++)
    {
      for(int j=0; j<numberOfQuarks; j++)
	{
	  if(Target.nucleonList.at(i).quarkList.at(j).collided == 1)
	    positionsQTW << Target.nucleonList.at(i).x + Target.nucleonList.at(i).quarkList.at(j).x << " " 
			 << Target.nucleonList.at(i).y + Target.nucleonList.at(i).quarkList.at(j).y << " " 
			 << Target.nucleonList.at(i).z + Target.nucleonList.at(i).quarkList.at(j).z << endl;
	  if(Target.nucleonList.at(i).quarkList.at(j).collided == 0)
	    positionsQTS << Target.nucleonList.at(i).x + Target.nucleonList.at(i).quarkList.at(j).x << " " 
			 << Target.nucleonList.at(i).y + Target.nucleonList.at(i).quarkList.at(j).y << " " 
			 << Target.nucleonList.at(i).z + Target.nucleonList.at(i).quarkList.at(j).z << endl;
	}
    }

 for(int i=0; i<Projectile.A; i++)
    {
      for(int j=0; j<numberOfQuarks; j++)
	{
	  if(Projectile.nucleonList.at(i).quarkList.at(j).collided == 1)
	    positionsQPW << Projectile.nucleonList.at(i).x + Projectile.nucleonList.at(i).quarkList.at(j).x << " " 
			 << Projectile.nucleonList.at(i).y + Projectile.nucleonList.at(i).quarkList.at(j).y << " " 
			 << Projectile.nucleonList.at(i).z + Projectile.nucleonList.at(i).quarkList.at(j).z << endl;
	  if(Projectile.nucleonList.at(i).quarkList.at(j).collided == 0)
	    positionsQPS << Projectile.nucleonList.at(i).x + Projectile.nucleonList.at(i).quarkList.at(j).x << " " 
			 << Projectile.nucleonList.at(i).y + Projectile.nucleonList.at(i).quarkList.at(j).y << " " 
			 << Projectile.nucleonList.at(i).z + Projectile.nucleonList.at(i).quarkList.at(j).z << endl;
	}
    }

  positionsQTW.close();
  positionsQTS.close();
  positionsQPW.close();
  positionsQPS.close();

}



void Glauber::computeSigmaNN(Random* random)
{
  int runs = 1000;
  int events = 5000;
  int counter = 0;
  double mean = 0.;
  double meanSq = 0.;
  double xShift, yShift;
  
  for(int j=0; j<runs; j++)
    {
      if (j%10==0)
	cout << j << endl;
      counter = 0;
      for(int i=0; i<events; i++)
	{
	  makeNuclei(random);
	  Projectile.nucleonList.at(0).collided = 0;
	  xShift = -5. + 10.*random->genrand64_real1();
	  yShift = -5. + 10.*random->genrand64_real1();
	  Projectile.nucleonList.at(0).x = xShift;
	  Projectile.nucleonList.at(0).y = yShift;
	  collide(random);
	  
	  if (Projectile.nucleonList.at(0).collided==1)
	    counter++;
	}
      mean += double(counter)/double(events) * 100. *10;
      meanSq += double(counter)/double(events) * 100. *10 * double(counter)/double(events) * 100. *10;
    }
  
  mean/=double(runs);
  meanSq/=double(runs);
  
  cout << mean << " mb +/- " << sqrt((meanSq-mean*mean))/sqrt(double(runs)) << " mb" << endl;
}
