#include "nrqcd.h"

namespace nrqcd{


double singlet(const double p, const double phip, const double k1, const double phik1,
               const double kprime, const double phikprime, const double k, const double phik, 
               const double m){
    
    double cosk1k = cos(phik1 - phik);
    double cospk = cos(phip - phik);
    double cospk1 = cos(phip - phik1);
    double cospkprime = cos(phip - phikprime);
    double cosk1kprime = cos(phik1 - phikprime);
    double kappa2 = k1*k1 + 4.*m*m;
    // l_T = k_T - p_T/2 + k_1T /2
    double lt2 = k*k + p*p/4. + k1*k1/4. - k*p*cospk + k*k1*cosk1k - p*k1*cospk1/2.;
    // l'_T = k'_T - p_T/2 + k_1T/2
    double lprimet2 = kprime*kprime + p*p/4. + k1*k1/4. - kprime*p*cospkprime + kprime*k1*cosk1kprime- p*k1*cospk1/2.;

    return k1*k1*kappa2*pow((1./(lt2+kappa2/4.)-1./(lprimet2+kappa2/4.)),2.)/(12.*m);

}

double octets10(const double p, const double phip, const double k1, const double phik1,
               const double k, const double phik, const double m){

    double cosk1k = cos(phik1 - phik);
    double cospk = cos(phip - phik);
    double cospk1 = cos(phip - phik1);
    double kappa2 = k1*k1 + 4*m*m;
    // l_T = k_T - p_T/2 + k_1T /2
    double lt2 = k*k + p*p/4. + k1*k1/4. - k*p*cospk + k*k1*cosk1k - p*k1*cospk1/2.;
    double kdotl = k*k - k*p*cospk/2. + k*k1*cosk1k/2.;
    double pdotk1 = p*k1*cospk1;
    double k1dotl = k1*k*cosk1k - k1*p*cospk1/2. + k1*k1/2.;
    double pdotl = p*k*cospk - p*p/2. + p*k1*cospk1/2.;

    return  2.*(k1*k1*lt2-k1dotl*k1dotl)/(m*pow(lt2+kappa2/4.,2.)); // changed to k1 dot l according to (B18f) in https://arxiv.org/pdf/1309.7337.pdf

 
 } 

 double octets13(const double p, const double phip, const double k1, const double phik1,
               const double k, const double phik, const double m){

    double cosk1k = cos(phik1 - phik);
    double cospk = cos(phip - phik);
    double cospk1 = cos(phip - phik1);
    double kappa2 = k1*k1 + 4*m*m;
    // l_T = k_T - p_T/2 + k_1T /2
    double lt2 = k*k + p*p/4. + k1*k1/4. - k*p*cospk + k*k1*cosk1k - p*k1*cospk1/2.;
    double kdotl = k*k - k*p*cospk/2. + k*k1*cosk1k/2.;
    double pdotk1 = p*k1*cospk1;
    double k1dotl = k1*k*cosk1k - k1*p*cospk1/2. + k1*k1/2.;
    double pdotl = p*k*cospk - p*p/2. + p*k1*cospk1/2.;


    double octet_two_a = 2.*k1*k1*(p*p+k1*k1-2.*pdotk1+4.*m*m)/(3.*m*m*m*(p*p+4.*m*m));
    double octet_two_b = -4.*k1*k1*(p*p+k1*k1-2.*pdotk1+pdotk1+4.*m*m)/(3.*m*(lt2+kappa2/4.)*(p*p+4.*m*m));
    double octet_two_c = k1*k1*kappa2/(6.*m*pow(lt2+kappa2/4.,2.));

    return octet_two_a+octet_two_b+octet_two_c;
 
 } 

  double octetp3j(const double p, const double phip, const double k1, const double phik1,
               const double k, const double phik, const double m){

    double cosk1k = cos(phik1 - phik);
    double cospk = cos(phip - phik);
    double cospk1 = cos(phip - phik1);
    double kappa2 = k1*k1 + 4*m*m;
    // l_T = k_T - p_T/2 + k_1T /2
    double lt2 = k*k + p*p/4. + k1*k1/4. - k*p*cospk + k*k1*cosk1k - p*k1*cospk1/2.;
    double kdotl = k*k - k*p*cospk/2. + k*k1*cosk1k/2.;
    double pdotk1 = p*k1*cospk1;
    double k1dotl = k1*k*cosk1k - k1*p*cospk1/2. + k1*k1/2.;
    double pdotl = p*k*cospk - p*p/2. + p*k1*cospk1/2.;
    
    double octet_three_a = (4.*k1*k1*lt2-2.*k1dotl*k1dotl)/(9.*m*m*m*pow(lt2+kappa2/4.,2.));
    double octet_three_b = (2.*k1*k1*k1dotl*(pdotl-k1dotl))/(9.*m*m*m*pow(lt2+kappa2/4.,3.));
    double octet_three_c = -(8.*m*m*(k1*k1*lt2-k1dotl*k1dotl))/(9.*m*m*m*pow(lt2+kappa2/4.,3.));
    double octet_three_d = k1*k1*kappa2*(pow(pdotl-k1dotl,2.)+4.*m*m*lt2)/(18.*m*m*m*pow(lt2+kappa2/4.,4.));
    
    return octet_three_a+octet_three_b+octet_three_c+octet_three_d;

 
 } 

}
