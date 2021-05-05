#ifndef NRQCD_H_
#define NRQCD_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>
#include <memory>
#include <iostream>

namespace nrqcd{
    
    double singlet(const double p, const double phip, const double k1, const double phik1,
                   const double kprime, const double phikprime, const double k, const double phik, 
                   const double m);
    
    double octets10(const double p, const double phip, const double k1, const double phik1,
               const double k, const double phik, const double m);

    double octets13(const double p, const double phip, const double k1, const double phik1,
               const double k, const double phik, const double m);

    double octetp3j(const double p, const double phip, const double k1, const double phik1,
               const double k, const double phik, const double m);              

}


#endif // NRQCD_H_