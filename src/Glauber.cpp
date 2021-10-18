#include "Glauber.h"

using namespace std;

Glauber::Glauber(Parameters *inParam, double inwidth)
{
  //  param = new Parameters();
  param = inParam; 
  gslRan = gsl_rng_alloc (gsl_rng_taus);
  numberOfQuarks = 3;
  width = inwidth;
}

// destructor
Glauber::~Glauber()
{
  //delete param;
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
	  if(theta<PI/40.*(j+1) && theta>PI/40.*j)
	    thetabins[j]+=1./samples/(PI/40.);
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
      thetahist << PI/40.*(j+0.5) << " " << thetabins[j] << endl;
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

  //  int rank = MPI::COMM_WORLD.Get_rank(); //number of current processor
} 

int Glauber::getA()
{
  return Target.A;
}

void Glauber::makeNuclei(Random *random, double Bp, double Bq, double Bqwidth, int Nq, double Yg, double YJPsi1, double YJPsi2, int Ydep)
{
  int Nx = param->getOutputNumberOfTransverseCells();

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
      
    }
  // if( Projectile.A == 1 )
  //   {
  //     Nucleon nucleon;
  //     nucleon.x = 0.;
  //     nucleon.y = 0.;
  //     nucleon.z = 0.;
  //     nucleon.collided = 0;
  //     Projectile.nucleonList.push_back(nucleon);
  //   }
  // else
  //   {
  //     for(int i=0; i<Projectile.A; i++)
  //       {
  //         Projectile.nucleonList.push_back(sampleRho(&Projectile, random));
  //       }
      
  //     //shift to center of mass
  //     double meanx=0, meany=0, meanz=0;
      
  //     for(int i=0; i<Projectile.A; i++)
  //       {
  //         meanx += Projectile.nucleonList.at(i).x;
  //         meany += Projectile.nucleonList.at(i).y;
  //         meanz += Projectile.nucleonList.at(i).z;
  //       }
      
  //     meanx /= double(Projectile.A);
  //     meany /= double(Projectile.A);
  //     meanz /= double(Projectile.A);
      
  //     for(int i=0; i<Projectile.A; i++)
  //       {
  //         Projectile.nucleonList.at(i).x -= meanx;
  //         Projectile.nucleonList.at(i).y -= meany;
  //         Projectile.nucleonList.at(i).z -= meanz;
  //       }     
  //   }



  // // if using quarks determine the quark positions within each nucleon
  // if(param->getUseQuarks())
  //   {
  //     double xQuark[numberOfQuarks];
  //     if( Target.A == 1 )
  //       {
  //         for(int j=0; j<numberOfQuarks; j++)
  //           {
  //             Target.nucleonList.at(0).quarkList.push_back(sampleQuark(random));
  //           }
  //       }
  //     else
  //       {
  //         for(int i=0; i<Target.A; i++)
  //           {
  //             for(int j=0; j<numberOfQuarks; j++)
  //       	{
  //       	  Target.nucleonList.at(i).quarkList.push_back(sampleQuark(random));
  //       	}
  //           }
  //       } 
  
  generateNucleusTA(&Target, random, Bp, Bq, Bqwidth,Nq, Yg, YJPsi1, YJPsi2, Ydep); 
  generateProtonTp(&Target, random, Bp, Bq, Bqwidth,Nq, Yg, YJPsi1, YJPsi2, Ydep); 
}

void Glauber::generateProtonTp(Nucleus *nuc, Random *random, double Bp, double Bq, double Bqwidth, int Nq, double Yg, double YJPsi1, double YJPsi2, int Ydep){
  // Bp, Bq are in GeV^-2
  double hbarc = 0.1973269804;
  
  double BqGauss;
  
  double gauss[Nq];
  double xq[Nq];
  double yq[Nq];

  for (int i = 0; i < Nq; i++) {
    gauss[i] = (exp(random->Gauss(0, width))) /
      std::exp(width * width / 2.0);
  }

  // placed here makes one size fluctuation for all hot spots simultaneously
  BqGauss= (exp(random->Gauss(0., Bqwidth))) /
    std::exp(Bqwidth * Bqwidth / 2.0);

  double YBp;
  //double avgxq = 0.;
  //double avgyq = 0.;
  // for (int iq = 0; iq < Nq; iq++) {
  //   avgxq += xq[iq];
  //   avgyq += yq[iq];
  // }
  // for (int iq = 0; iq < Nq; iq++) {
  //   xq[iq] -= avgxq / double(Nq);
  //   yq[iq] -= avgyq / double(Nq);
  // }
 
  if(Ydep==0){ 
    for (int iq = 0; iq < Nq; iq++) {
      xq[iq] = sqrt(Bp * hbarc * hbarc) * random->Gauss();
      yq[iq] = sqrt(Bp * hbarc * hbarc) * random->Gauss();
    }
    for(int ix=0; ix<40; ix++){
      double x = (double(ix)/40.*4.-2.);
      for(int iy=0; iy<40; iy++){
        double y = (double(iy)/40.*4.-2.);
        Tpgrid2D[0][ix][iy] = 0.;
        Tpgrid2D[1][ix][iy] = 0.;
        Tpgrid2D[2][ix][iy] = 0.;
        for(int i=0; i<Nq; i++){
         Tpgrid2D[0][ix][iy] += exp(-((x-xq[i])*(x-xq[i])+(y-yq[i])*(y-yq[i]))/hbarc/hbarc/2./(Bq*BqGauss))*gauss[i]/(Bq*BqGauss)/double(Nq)*((Bq*BqGauss)+Bp); // normalize to 1 at zero on average
          Tpgrid2D[1][ix][iy] += exp(-((x-xq[i])*(x-xq[i])+(y-yq[i])*(y-yq[i]))/hbarc/hbarc/2./(Bq*BqGauss))*gauss[i]/(Bq*BqGauss)/double(Nq)*((Bq*BqGauss)+Bp); // normalize to 1 at zero on average
          Tpgrid2D[2][ix][iy] += exp(-((x-xq[i])*(x-xq[i])+(y-yq[i])*(y-yq[i]))/hbarc/hbarc/2./(Bq*BqGauss))*gauss[i]/(Bq*BqGauss)/double(Nq)*((Bq*BqGauss)+Bp); // normalize to 1 at zero on average
        }
      }
    }
  }
  else{
    YBp=Bp*(0.15 + 0.042*pow((Yg - 4.6),2.));
    for (int iq = 0; iq < Nq; iq++) {
      xq[iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
      yq[iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
    }
    double YBq = Bq*(0.15 + 0.042*pow((Yg - 4.6),2.));
    for(int ix=0; ix<40; ix++){
      double x = (double(ix)/40.*4.-2.);
      for(int iy=0; iy<40; iy++){
        double y = (double(iy)/40.*4.-2.);
        Tpgrid2D[0][ix][iy] = 0.;
        for(int i=0; i<Nq; i++){
          Tpgrid2D[0][ix][iy] += exp(-((x-xq[i])*(x-xq[i])+(y-yq[i])*(y-yq[i]))/hbarc/hbarc/2./(YBq*BqGauss))*gauss[i]/(YBq*BqGauss)/double(Nq)*((YBq*BqGauss)+YBp); // normalize to 1 at zero on average
        }
      }
    }
    YBp=Bp*(0.15 + 0.042*pow((YJPsi1 - 4.6),2.));
    for (int iq = 0; iq < Nq; iq++) {
      xq[iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
      yq[iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
    }
    YBq = Bq*(0.15 + 0.042*pow((YJPsi1 - 4.6),2.));
    for(int ix=0; ix<40; ix++){
      double x = (double(ix)/40.*4.-2.);
      for(int iy=0; iy<40; iy++){
        double y = (double(iy)/40.*4.-2.);
        Tpgrid2D[1][ix][iy] = 0.;
        for(int i=0; i<Nq; i++){
          Tpgrid2D[1][ix][iy] += exp(-((x-xq[i])*(x-xq[i])+(y-yq[i])*(y-yq[i]))/hbarc/hbarc/2./(YBq*BqGauss))*gauss[i]/(YBq*BqGauss)/double(Nq)*((YBq*BqGauss)+YBp); // normalize to 1 at zero on average
        }
      }
    }
    YBp=Bp*(0.15 + 0.042*pow((YJPsi2 - 4.6),2.));
    for (int iq = 0; iq < Nq; iq++) {
      xq[iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
      yq[iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
    }
    YBq=Bq*(0.15 + 0.042*pow((YJPsi2 - 4.6),2.));
    for(int ix=0; ix<40; ix++){
      double x = (double(ix)/40.*4.-2.);
      for(int iy=0; iy<40; iy++){
        double y = (double(iy)/40.*4.-2.);
        Tpgrid2D[2][ix][iy] = 0.;
        for(int i=0; i<Nq; i++){
          Tpgrid2D[2][ix][iy] += exp(-((x-xq[i])*(x-xq[i])+(y-yq[i])*(y-yq[i]))/hbarc/hbarc/2./(YBq*BqGauss))*gauss[i]/(YBq*BqGauss)/double(Nq)*((YBq*BqGauss)+YBp); // normalize to 1 at zero on average
        }
      }
    }
  }
}

void Glauber::generateNucleusTA(Nucleus *nuc, Random *random, double Bp, double Bq, double Bqwidth, int Nq, double Yg, double YJPsi1, double YJPsi2, int Ydep){
  // Bp is in GeV^-2
  double hbarc = 0.1973269804;
  int useQuarks = 1;
  double BqGauss;
  // stringstream strfilename;
  // strfilename << "TA.dat";
  // string filename;
  // filename = strfilename.str();
  // fstream fout(filename.c_str(), ios::out);
  
  double gauss[nuc->nucleonList.size()][Nq];
  double xq[nuc->nucleonList.size()][Nq];
  double yq[nuc->nucleonList.size()][Nq];

  for (unsigned int i = 0; i < nuc->nucleonList.size(); i++) {
    for (int iq = 0; iq < Nq; iq++) {
      gauss[i][iq] = (exp(random->Gauss(0, width))) /
        std::exp(width * width / 2.0);
    }
  }

  // placed here makes one size fluctuation for all hot spots simultaneously
  BqGauss= (exp(random->Gauss(0., Bqwidth))) /
    std::exp(Bqwidth * Bqwidth / 2.0);
  
  double YBp;

  if (useQuarks == 1){
    for (unsigned int i = 0; i < nuc->nucleonList.size(); i++) {
      //      double avgxq = 0.;
      //double avgyq = 0.;
      for (int iq = 0; iq < Nq; iq++) {
        xq[i][iq] = sqrt(Bp * hbarc * hbarc) * random->Gauss();
        yq[i][iq] = sqrt(Bp * hbarc * hbarc) * random->Gauss();
      }
      // for (int iq = 0; iq < Nq; iq++) {
      //   avgxq += xq[i][iq];
      //   avgyq += yq[i][iq];
      // }
      // for (int iq = 0; iq < Nq; iq++) {
      //   xq[i][iq] -= avgxq / double(Nq);
      //   yq[i][iq] -= avgyq / double(Nq);
      // }
    }
  }
  // //  double gauss[nuc->nucleonList.size()];

  // // width=0.5 for now
  // for (unsigned int i = 0; i < nuc->nucleonList.size(); i++) {
  //   //   gauss[i] = (exp(random->Gauss(0, 0.5))) /
  //   //  std::exp(0.5 * 0.5 / 2.0);
  //   gauss[i] = (exp(random->Gauss(0, 0.5))) /
  //     std::exp(0.5 * 0.5 / 2.0);
  // }

  if(Ydep==0){
    for(int ix=0; ix<200; ix++){
      double x = (double(ix)/200.*20.-10.);
      for(int iy=0; iy<200; iy++){
        double y = (double(iy)/200.*20.-10.);
        TAgrid2D[0][ix][iy] = 0.;
        TAgrid2D[1][ix][iy] = 0.;
        TAgrid2D[2][ix][iy] = 0.;
        for(unsigned int i=0; i<nuc->nucleonList.size(); i++){
          if(useQuarks == 1){
            for (int iq = 0; iq < Nq; iq++) {
              double xpos = nuc->nucleonList.at(i).x+xq[i][iq]; 
              double ypos = nuc->nucleonList.at(i).y+yq[i][iq]; 
              TAgrid2D[0][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))
                                         /hbarc/hbarc/2./(Bq*BqGauss))*gauss[i][iq]/(Bq*BqGauss)/double(Nq)*((Bq*BqGauss)+Bp); // normalize proton to 1 at zero on average
              TAgrid2D[1][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))
                                         /hbarc/hbarc/2./(Bq*BqGauss))*gauss[i][iq]/(Bq*BqGauss)/double(Nq)*((Bq*BqGauss)+Bp); // normalize proton to 1 at zero on average
              TAgrid2D[2][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))
                                         /hbarc/hbarc/2./(Bq*BqGauss))*gauss[i][iq]/(Bq*BqGauss)/double(Nq)*((Bq*BqGauss)+Bp); // normalize proton to 1 at zero on average
            }
          }
          else{
            double xpos = nuc->nucleonList.at(i).x; 
            double ypos = nuc->nucleonList.at(i).y; 
            TAgrid2D[0][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))/hbarc/hbarc/2./Bp)*gauss[i][0];
            TAgrid2D[1][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))/hbarc/hbarc/2./Bp)*gauss[i][0];
            TAgrid2D[2][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))/hbarc/hbarc/2./Bp)*gauss[i][0];
          }
        }
        //fout << x << " " << y << " " << TAgrid2D[ix][iy] << endl;
      }
    }
  }
  else{
    YBp=Bp*(0.15 + 0.042*pow((-Yg - 4.6),2.));
    for (unsigned int i = 0; i < nuc->nucleonList.size(); i++) {
      for (int iq = 0; iq < Nq; iq++) {
        xq[i][iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
        yq[i][iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
      }
    }
    double YBq = Bq*(0.15 + 0.042*pow((-Yg - 4.6),2.));
    for(int ix=0; ix<200; ix++){
      double x = (double(ix)/200.*20.-10.);
      for(int iy=0; iy<200; iy++){
        double y = (double(iy)/200.*20.-10.);
        TAgrid2D[0][ix][iy] = 0.;
        for(unsigned int i=0; i<nuc->nucleonList.size(); i++){
          if(useQuarks == 1){
            for (int iq = 0; iq < Nq; iq++) {
              double xpos = nuc->nucleonList.at(i).x+xq[i][iq]; 
              double ypos = nuc->nucleonList.at(i).y+yq[i][iq]; 
              TAgrid2D[0][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))
                                         /hbarc/hbarc/2./(YBq*BqGauss))*gauss[i][iq]/(YBq*BqGauss)/double(Nq)*((YBq*BqGauss)+YBp); // normalize proton to 1 at zero on average
            }
          }
          else{
            double xpos = nuc->nucleonList.at(i).x; 
            double ypos = nuc->nucleonList.at(i).y; 
            TAgrid2D[0][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))/hbarc/hbarc/2./YBp)*gauss[i][0];
          }
        }
      }
    }
    YBp=Bp*(0.15 + 0.042*pow((-YJPsi1 - 4.6),2.));
    for (unsigned int i = 0; i < nuc->nucleonList.size(); i++) {
      for (int iq = 0; iq < Nq; iq++) {
        xq[i][iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
        yq[i][iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
      }
    }
    YBq = Bq*(0.15 + 0.042*pow((-YJPsi1 - 4.6),2.));
    for(int ix=0; ix<200; ix++){
      double x = (double(ix)/200.*20.-10.);
      for(int iy=0; iy<200; iy++){
        double y = (double(iy)/200.*20.-10.);
        TAgrid2D[1][ix][iy] = 0.;
        for(unsigned int i=0; i<nuc->nucleonList.size(); i++){
          if(useQuarks == 1){
            for (int iq = 0; iq < Nq; iq++) {
              double xpos = nuc->nucleonList.at(i).x+xq[i][iq]; 
              double ypos = nuc->nucleonList.at(i).y+yq[i][iq]; 
              TAgrid2D[1][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))
                                         /hbarc/hbarc/2./(YBq*BqGauss))*gauss[i][iq]/(YBq*BqGauss)/double(Nq)*((YBq*BqGauss)+YBp); // normalize proton to 1 at zero on average
            }
          }
          else{
            double xpos = nuc->nucleonList.at(i).x; 
            double ypos = nuc->nucleonList.at(i).y; 
            TAgrid2D[1][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))/hbarc/hbarc/2./YBp)*gauss[i][0];
          }
        }
      }
    }
    YBp=Bp*(0.15 + 0.042*pow((-YJPsi2 - 4.6),2.));
    for (unsigned int i = 0; i < nuc->nucleonList.size(); i++) {
      for (int iq = 0; iq < Nq; iq++) {
        xq[i][iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
        yq[i][iq] = sqrt(YBp * hbarc * hbarc) * random->Gauss();
      }
    }
    YBq = Bq*(0.15 + 0.042*pow((-YJPsi2 - 4.6),2.));
    for(int ix=0; ix<200; ix++){
      double x = (double(ix)/200.*20.-10.);
      for(int iy=0; iy<200; iy++){
        double y = (double(iy)/200.*20.-10.);
        TAgrid2D[2][ix][iy] = 0.;
        for(unsigned int i=0; i<nuc->nucleonList.size(); i++){
          if(useQuarks == 1){
            for (int iq = 0; iq < Nq; iq++) {
              double xpos = nuc->nucleonList.at(i).x+xq[i][iq]; 
              double ypos = nuc->nucleonList.at(i).y+yq[i][iq]; 
              TAgrid2D[2][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))
                                         /hbarc/hbarc/2./(YBq*BqGauss))*gauss[i][iq]/(YBq*BqGauss)/double(Nq)*((YBq*BqGauss)+YBp); // normalize proton to 1 at zero on average
            }
          }
          else{
            double xpos = nuc->nucleonList.at(i).x; 
            double ypos = nuc->nucleonList.at(i).y; 
            TAgrid2D[2][ix][iy] += exp(-((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos))/hbarc/hbarc/2./YBp)*gauss[i][0];
          }
        }
      }
    }
  }
}  

//takes x and y in GeV^-1
 double Glauber::returnNucleusTA(double x, double y, int Ybin){
  double hbarc = 0.1973269804;
  // simple linear interpolation
  
  if (x*hbarc<-9.9 || x*hbarc>9.9){
    return 0.;
  }      
  if (y*hbarc<-9.9 || y*hbarc>9.9){
    return 0.;
  }      
 
  int ix = int((x*hbarc+10.)*200./20.);
  int iy = int((y*hbarc+10.)*200./20.);

  if (ix<0 || ix>198){
    cout << " x out of range: " << endl;
    cout << "x = " << x*hbarc << endl;
    cout << "ix = " << ix << endl;
  }
  if (iy<0 || iy>198){
    cout << " y out of range: " << endl;
    cout << "y = " << y*hbarc << endl;
    cout << "iy = " << iy << endl;
    return 0.;
  }

  double TAx1 = (double(ix+1)-(x*hbarc+10.)*200./20.)*TAgrid2D[Ybin][ix][iy] + ((x*hbarc+10.)*200./20.-double(ix))*TAgrid2D[Ybin][ix+1][iy];
  double TAx2 = (double(ix+1)-(x*hbarc+10.)*200./20.)*TAgrid2D[Ybin][ix][iy+1] + ((x*hbarc+10.)*200./20.-double(ix))*TAgrid2D[Ybin][ix+1][iy+1];
  double TA = (double(iy+1)-(y*hbarc+10.)*200./20.)*TAx1 + ((y*hbarc+10.)*200./20.-double(iy))*TAx2;

  return TA;
}

//takes x and y in GeV^-1
 double Glauber::returnProtonTp(double x, double y, int Ybin){
  double hbarc = 0.1973269804;
  // simple linear interpolation
  
  if (x*hbarc<-1.9 || x*hbarc>1.9){
    return 0.;
  }      
  if (y*hbarc<-1.9 || y*hbarc>1.9){
    return 0.;
  }      
 
  int ix = int((x*hbarc+2.)*40./4.);
  int iy = int((y*hbarc+2.)*40./4.);

  if (ix<0 || ix>38){
    cout << " x out of range: " << endl;
    cout << "x = " << x*hbarc << endl;
    cout << "ix = " << ix << endl;
  }
  if (iy<0 || iy>38){
    cout << " y out of range: " << endl;
    cout << "y = " << y*hbarc << endl;
    cout << "iy = " << iy << endl;
    return 0.;
  }

  double Tpx1 = (double(ix+1)-(x*hbarc+2.)*40./4.)*Tpgrid2D[Ybin][ix][iy] + ((x*hbarc+2.)*40./4.-double(ix))*Tpgrid2D[Ybin][ix+1][iy];
  double Tpx2 = (double(ix+1)-(x*hbarc+2.)*40./4.)*Tpgrid2D[Ybin][ix][iy+1] + ((x*hbarc+2.)*40./4.-double(ix))*Tpgrid2D[Ybin][ix+1][iy+1];
  double Tp = (double(iy+1)-(y*hbarc+2.)*40./4.)*Tpx1 + ((y*hbarc+2.)*40./4.-double(iy))*Tpx2;

  return Tp;
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

Nucleon Glauber::sampleRho(Nucleus *nuc, Random *random)
{
  Nucleon nucleon;
  double r, x, y, z, tmp;
  double phi, theta;
  double a_WS = nuc->a_WS;
  double R_WS = nuc->R_WS;
   
  phi = 2.*PI*random->genrand64_real1();
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
  //z=r*cos(theta);
  
  // set values of nucleon
  nucleon.x=x;
  nucleon.y=y;
  nucleon.z=0.; // z
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

  phi = 2.*PI*random->genrand64_real1();
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
  //  z=r*cos(theta);
  
  // set values of nucleon
  quark.x=x;
  quark.y=y;
  quark.z=0.;
  quark.collided=0;

  return quark; 
}

