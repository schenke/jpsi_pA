#include "Parameters.h"

string Parameters::strWord(int index, string line)
{
  int count = 0; // number of read words
  string word; // the resulting word
  for (long unsigned int i = 0 ; i < line.length(); i++) 
    { // iterate over all characters in 'line'
      if (line[i] == ' ')
	{ // if this character is a space we might be done reading a word from 'line'
	  if (line[i+1] != ' ') 
	    { // next character is not a space, so we are done reading a word
	      count++; // increase number of read words
	      if (count == index) 
		{ // was this the word we were looking for?
		  return word;
		}
	      word.clear();
	      word = ""; // nope it wasn't .. so reset word and start over with the next one in 'line'
	    }
	}
      else 
	{ // not a space .. so append the character to 'word'
	  word += line[i];
	}
    }
  return "0";
}

// read in parameter from the input file
bool Parameters::readParameter(string parameter)
{
  bool success=0;
  ifstream inputFile("input");
  string line, variable, value;
  while ( getline (inputFile,line) )
   {
     line += ' ';
     if (strWord(1, line) == parameter )
       { 
	 success = 1;
	 temp = strWord(2, line);
       }
   }
  inputFile.close();
  return success;
}

// set parameters using defaults and the input file
void Parameters::setParameters()
{
  vector<string> parameters;
  stringstream convert;
  double doubleParameter;
  int intParameter;
  long longParameter;
 
  // set default first          
  // read from file if present
  
  setProjectile("p");     
  if(readParameter("Projectile")) setProjectile(temp);
 
  setTarget("Pb");         
  if(readParameter("Target")) setTarget(temp);
  
  setb(0.);                
  if(readParameter("b")){ convert<<temp; convert >> doubleParameter; setb(doubleParameter); 
    convert.str(""); convert.clear(); } 

  setTimeForSeed(1);  
  if(readParameter("timeForSeed")){ convert<<temp; convert >> intParameter; setTimeForSeed(intParameter); 
    convert.str(""); convert.clear(); } 

 
  setSeed(1);              
  if(readParameter("seed"))
    {
      if (getTimeForSeed())
	{ convert<<temp; convert >> longParameter; setSeed(time(0)+longParameter*10000+1000); 
	  convert.str(""); convert.clear(); }
      else 
	{ convert<<temp; convert >> longParameter; setSeed(longParameter*10000+1000); 
	  convert.str(""); convert.clear(); }
    }

  
  setGaussianWounding(0);  
  if(readParameter("gaussianWounding")){ convert<<temp; convert >> intParameter; setGaussianWounding(intParameter); convert.str(""); 
    convert.clear(); }  
 
  setUseQuarks(0);  
  if(readParameter("useQuarks")){ convert<<temp; convert >> intParameter; setUseQuarks(intParameter); convert.str(""); 
    convert.clear(); } 


  cout << "1" << endl;



  
  setOutputNumberOfTransverseCells(200);  
  if(readParameter("outputNumberOfTransverseCells")){ convert<<temp; convert >> intParameter; setOutputNumberOfTransverseCells(intParameter); convert.str(""); 
    convert.clear(); } 

  setUseEnergyDependentCrossSection(0);  
  if(readParameter("useEnergyDependentCrossSection")){ convert<<temp; convert >> intParameter; setUseEnergyDependentCrossSection(intParameter); convert.str(""); 
    convert.clear(); } 
  
  setSigmaNN(42.);          
  if(readParameter("sigmaNN")){ convert<<temp; convert >> doubleParameter; setSigmaNN(doubleParameter); 
    convert.str(""); convert.clear(); } 
 
  setSigmaQQ(9.36);          
  if(readParameter("sigmaQQ")){ convert<<temp; convert >> doubleParameter; setSigmaQQ(doubleParameter); 
    convert.str(""); convert.clear(); } 

   setRoots(200.);          
  if(readParameter("roots")){ convert<<temp; convert >> doubleParameter; setRoots(doubleParameter); 
    convert.str(""); convert.clear(); } 

  
  // when using energy dependent cross sections use the fit values for the nucleon nucleon cross sections and use a function that is proportional
  // to that for the qq cross section with the normalization such that the entered sigmaQQ is the value at roots = 200 GeV 

  if(getUseEnergyDependentCrossSection() == 1)
    {
      //use result from http://arxiv.org/pdf/hep-ph/0206028.pdf


   
      //     setSigmaNNinel((25.2-0.05*log(getRoots())+0.56*log(getRoots())*log(getRoots()))+45.2/pow(getRoots(),0.9)+33.8/pow(getRoots(),1.1));
      //   setSigmaQQinel(getSigmaQQ()*((25.2-0.05*log(getRoots())+0.56*log(getRoots())*log(getRoots()))+45.2/pow(getRoots(),0.9)+33.8/pow(getRoots(),1.1)));
    

      setSigmaNN((44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots())));
      setSigmaNNinel((44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots()))
		     -(11.4-1.52*log(getRoots()*getRoots())+0.13*log(getRoots()*getRoots())*log(getRoots()*getRoots()))); // total - elastic

      setSigmaQQ(getSigmaQQ()*(44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots()))/50.725019604018556);
      setSigmaQQinel(getSigmaQQ()*((44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots()))
				   -(11.4-1.52*log(getRoots()*getRoots())+0.13*log(getRoots()*getRoots())*log(getRoots()*getRoots())))
       		     /(44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots())));
 
      // setSigmaNN((47-3.6*log(getRoots()*getRoots())+0.39*log(getRoots()*getRoots())*log(getRoots()*getRoots())));
      // setSigmaNNinel((47-3.6*log(getRoots()*getRoots())+0.39*log(getRoots()*getRoots())*log(getRoots()*getRoots()))
      // 		 -(11.4-1.52*log(getRoots()*getRoots())+0.13*log(getRoots()*getRoots())*log(getRoots()*getRoots()))); // total - elastic

      // setSigmaQQ(getSigmaQQ()*(47-3.6*log(getRoots()*getRoots())+0.39*log(getRoots()*getRoots())*log(getRoots()*getRoots()))/50.725019604018556);
      // setSigmaQQinel(getSigmaQQ()*((47-3.6*log(getRoots()*getRoots())+0.39*log(getRoots()*getRoots())*log(getRoots()*getRoots()))
      // 			       -(11.4-1.52*log(getRoots()*getRoots())+0.13*log(getRoots()*getRoots())*log(getRoots()*getRoots())))
      // 		     /(47-3.6*log(getRoots()*getRoots())+0.39*log(getRoots()*getRoots())*log(getRoots()*getRoots())));

      // adjust for stronger energy dependence in heavy ion collisions (see ALICE http://arxiv.org/pdf/1011.3916v3.pdf)
      //setSigmaNN(getSigmaNN()*(pow(getRoots()*getRoots(),0.05)/pow(200*200,0.05)));
      //setSigmaNNinel(getSigmaNNinel()*(pow(getRoots()*getRoots(),0.05)/pow(200*200,0.05)));
      //setSigmaQQ(getSigmaQQ()*(pow(getRoots()*getRoots(),0.05)/pow(200*200,0.05)));
      //setSigmaQQinel(getSigmaQQinel()*(pow(getRoots()*getRoots(),0.05)/pow(200*200,0.05)));
    }			     
  else
    {

   
      //     setSigmaNNinel((25.2-0.05*log(getRoots())+0.56*log(getRoots())*log(getRoots()))+45.2/pow(getRoots(),0.9)+33.8/pow(getRoots(),1.1));
      //setSigmaQQinel(getSigmaQQ()*((25.2-0.05*log(getRoots())+0.56*log(getRoots())*log(getRoots()))+45.2/pow(getRoots(),0.9)+33.8/pow(getRoots(),1.1)));

      setSigmaQQinel(getSigmaQQ()*((44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots()))
				   -(11.4-1.52*log(getRoots()*getRoots())+0.13*log(getRoots()*getRoots())*log(getRoots()*getRoots())))
       		     /(44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots()))); 
      setSigmaNNinel(getSigmaNN() *((44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots()))
       			       -(11.4-1.52*log(getRoots()*getRoots())+0.13*log(getRoots()*getRoots())*log(getRoots()*getRoots())))
       		     /(44.4-2.9*log(getRoots()*getRoots())+0.33*log(getRoots()*getRoots())*log(getRoots()*getRoots())));
      
      // setSigmaNNinel(getSigmaNNinel()*(pow(getRoots()*getRoots(),0.05)/pow(200*200,0.05)));
      // setSigmaQQinel(getSigmaQQinel()*(pow(getRoots()*getRoots(),0.05)/pow(200*200,0.05)));
    }

 setOutputMaximalTransverseSize(12.);
  if(readParameter("outputMaximalTransverseSize")){ convert<<temp; convert >> doubleParameter; setOutputMaximalTransverseSize(doubleParameter); 
    convert.str(""); convert.clear(); } 

 setTransverseGaussianSmearingWidth(0.4);          
  if(readParameter("transverseGaussianSmearingWidth")){ convert<<temp; convert >> doubleParameter; setTransverseGaussianSmearingWidth(doubleParameter); 
    convert.str(""); convert.clear(); } 


  int rank = 0;//MPI::COMM_WORLD.Get_rank(); //number of current processor
  if(rank==0)
    {
  // output parameter info
  cout << "------------------------------------------------------------------------------------------------" << endl;
  cout << "| Parameters: " << endl;
  cout << "------------------------------------------------------------------------------------------------" << endl;
  cout << setw(45) << left << "| Projectile: " << setw(20) << getProjectile() << setw(45) << " | (default p)" << endl;
  cout << setw(45) << left << "| Target: " << setw(20) << getTarget() << setw(45) << " | (default Pb)" << endl;
  convert << getb() << " fm";
  cout << setw(45) << left << "| b: " << setw(20) << convert.str() << setw(45) << " | (default 0 fm)" << endl; convert.str(""); convert.clear();
  cout << setw(45) << left << "| Gaussian Wounding: " << setw(20) << getGaussianWounding() << setw(45) << " | 0: off, 1: on (default 0)" << endl;
  cout << setw(45) << left << "| use quarks: " << setw(20) << getUseQuarks() << setw(45) << " | 0: no, 1: yes (default 0)" << endl;
  cout << setw(45) << left << "| use energy dependent cross section: " << setw(20) << getUseEnergyDependentCrossSection() << setw(45) << " | 0: no, 1: yes (default 0)" << endl;
  if(getUseQuarks())
    {
      convert << getSigmaQQ() << " mb";
      cout << setw(45) << left << "| sigmaQQ: " << setw(20) << convert.str() << setw(45) << " | (default 9.36 mb)" << endl; convert.str(""); convert.clear();
    } 
  else
    {
      convert << getSigmaNN() << " mb";
      cout << setw(45) << left << "| sigmaNN: " << setw(20) << convert.str() << setw(45) << " | (default 42 mb)" << endl; convert.str(""); convert.clear();
    }
  convert << getRoots() << " GeV";
  cout << setw(45) << left << "| center of mass energy: " << setw(20) << convert.str() << setw(45) << " | (default 200 GeV)" << endl; convert.str(""); convert.clear();
  cout << setw(45) << left << "| using time for seed: " << setw(20) << getTimeForSeed() << setw(45) << " | 0: no, 1: yes (default 1)" << endl;
  cout << setw(45) << left << "| seed: " << setw(20) << getSeed() << setw(45) << " |" << endl;
  cout << setw(45) << left << "| number of transverse cells in output grid: " << setw(20) << getOutputNumberOfTransverseCells() << setw(45) << " | (default 200)" << endl;
  convert << getOutputMaximalTransverseSize() << " fm";
  cout << setw(45) << left << "| maximal transverse size: " << setw(20) << convert.str() << setw(45) << " | (default 12.0 fm)" << endl; convert.str(""); convert.clear();
  convert << getTransverseGaussianSmearingWidth() << " fm";
  cout << setw(45) << left << "| Gaussian smearing width in x and y: " << setw(20) << convert.str() << setw(45) << " | (default 0.4 fm)" << endl; convert.str(""); convert.clear();
  cout << "------------------------------------------------------------------------------------------------" << endl;
  cout << endl;
    }    

}


