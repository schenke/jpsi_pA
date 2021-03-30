#ifndef HARD_H_
#define HARD_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>
#include <memory>

namespace Hard{
  
  double all(const double p, const double phip, const double q, const double phiq, 
              const double k1, const double phik1, const double k2, const double phik2, 
              const double k, const double phik, const double yp, const double yq, const double m);
  double qqqq(const double p, const double phip, const double q, const double phiq, 
              const double k1, const double phik1, const double k2, const double phik2, 
              const double k, const double phik, const double yp, const double yq, const double m);
  double qqqq_col(const double p, const double phip, const double q, const double phiq, 
              const double k1, const double phik1, const double k2, const double phik2, 
              const double k, const double phik, const double yp, const double yq, const double m);
  double qqg(const double p, const double phip, const double q, const double phiq, 
              const double k1, const double phik1, const double k2, const double phik2, 
              const double k, const double phik, const double yp, const double yq, const double m);
  double gg(const double p, const double phip, const double q, const double phiq, 
              const double k1, const double phik1, const double k2, const double phik2, 
              const double k, const double phik, const double yp, const double yq, const double m);

}

#endif // HARD_H_
