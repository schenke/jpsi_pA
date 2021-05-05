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

namespace NRQCD{
    double singlet(const double p, const double phip, const double k1, const double phik1,
                   const double kprime, const double phikprime, const double k, const double phik, 
                   const double m);
    
    double octet_s10(const double p, const double phip, const double k1, const double phik1,
                 const double kprime, const double phikprime, const double k, const double phik, 
                 const double m);

    double octet_s13(const double p, const double phip, const double k1, const double phik1,
                 const double kprime, const double phikprime, const double k, const double phik, 
                 const double m);

    double octet_p3j(const double p, const double phip, const double k1, const double phik1,
                 const double kprime, const double phikprime, const double k, const double phik, 
                 const double m);              

}


#endif // NRQCD_H_