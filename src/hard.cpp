#include "hard.h"

namespace Hard{
  
  // All momenta passed are transverse momenta. + and - components are computed from them and the rapidities, 
  // and m
  // I chose to pass all parameters to all functions, even if redundant, to avoid mistakes later.
  double all(const double p, const double phip, const double q, const double phiq, 
             const double k1, const double phik1, const double k2, const double phik2, 
             const double k, const double phik, const double yp, const double yq, const double m) {
    
    double cosqk = cos(phiq-phik);
    double cosqk1 = cos(phiq-phik1);
    double coskk1 = cos(phik-phik1);
    double cosk1p = cos(phik1-phip);
    double cosk1k2 = cos(phik1-phik2);
    double cosk2q = cos(phik2-phiq);
    
    double pplus = sqrt(m*m+p*p)*exp(yp)/sqrt(2.);
    double pminus = sqrt(m*m+p*p)*exp(-yp)/sqrt(2.);
    double qplus = sqrt(m*m+q*q)*exp(yq)/sqrt(2.);
    double qminus = sqrt(m*m+q*q)*exp(-yq)/sqrt(2.);

    // a_T = q_T-k_T (2D vectors)
    double at2 = q*q + k*k - 2*q*k*cosqk;
    // b_T = q_T-k_T-k1_T (2D vectors)
    double bt2 = at2 + k1*k1 -2.*q*k1*cosqk1 + 2.*k*k1*coskk1;
    double adotb = at2 + k*k1*coskk1 - q*k1*cosqk1;
    double adotq = q*q - k*q*cosqk;
    double bdotq = adotq - k1*q*cosqk1;
    double adotp = q*p*cos(phiq-phip) - k*p*cos(phik-phip);
    double bdotp = adotp - k1*p*cosk1p;
    
    // Lipatov vertex
    double Cplus = pplus + qplus - k1*k1/(pminus+qminus);
    double Cminus = k2*k2/(pplus+qplus) - pminus - qminus;
    double adotC = q*k2*cosk2q -q*k1*cosqk1
      -k*k2*cos(phik-phik2) +k*k1*coskk1;
    double bdotC = adotC - k1*k2*cosk1k2 + k1*k1;
    // four vector products
    double Cp = (Cplus*pminus + Cminus*pplus) 
      - k2*p*cos(phik2-phip) + k1*p*cosk1p;
    double Cq = (Cplus*qminus + Cminus*qplus) 
      - k2*q*cosk2q + k1*q*cosqk1;
    double pq = (pplus*qminus + pminus*qplus) 
      - p*q*cos(phip-phiq);
    double C2 = 2.*Cplus*Cminus - ( k2*k2+k1*k1-2.*k1*k2*cosk1k2 );
    return (32.*pplus*qplus*(m*m+at2)*(m*m+bt2)
            /pow((2.*pplus*(m*m+at2)+2.*qplus*(m*m+bt2)),2.))
      + (4./((m*m+pq)*(pplus*(m*m+at2)+qplus*(m*m+bt2)))*
         ((m*m+adotb)*(qplus*Cp+pplus*Cq-Cplus*(m*m+pq))
          +Cplus*((m*m+bdotq)*(m*m-adotp)-(m*m+adotq)*(m*m-bdotp))
          +pplus*(adotC*(m*m+bdotq)-bdotC*(m*m+adotq))
          +qplus*(adotC*(m*m-bdotp)-bdotC*(m*m-adotp))))
      +((2.*Cp*Cq - (m*m+pq)*C2)/pow((m*m+pq),2.));
  }

  double qqqq(const double p, const double phip, const double q, const double phiq, 
              const double k1, const double phik1, const double k2, const double phik2, 
              const double k, const double phik, const double yp, const double yq, const double m) {
    
    double pplus = sqrt(m*m+p*p)*exp(yp)/sqrt(2.);
    double pminus = sqrt(m*m+p*p)*exp(-yp)/sqrt(2.);
    double qplus = sqrt(m*m+q*q)*exp(yq)/sqrt(2.);
    double qminus = sqrt(m*m+q*q)*exp(-yq)/sqrt(2.);

    // a_T = q_T-k_T (2D vectors)
    double at2 = q*q + k*k - 2*q*k*cos(phiq-phik);
    // b_T = q_T-k_T-k1_T (2D vectors)
    double bt2 = at2 + k1*k1 -2.*q*k1*cos(phiq-phik1) + 2.*k*k1*cos(phik-phik1);

    return 32.*pplus*qplus*(m*m+at2)*(m*m+bt2)/pow((2.*pplus*(m*m+at2)+2.*qplus*(m*m+bt2)),2.);

  }

  // Collinear version of the hard matrix element, i.e. expansion at linear order in k1
  double qqqq_col(const double p, const double phip, const double q, const double phiq, 
            const double k1, const double phik1, const double k2, const double phik2, 
            const double k, const double phik, const double yp, const double yq, const double m){ 
  // a_T = q_T-k_T (2D vectors)
    double at2 = q*q + k*k - 2*q*k*cos(phiq-phik);
    double pplus = sqrt(m*m+p*p)*exp(yp);
    double qplus = sqrt(m*m+q*q)*exp(yq);

    double prefactor = 8.*pplus*qplus/(pow(pplus+qplus,2)*(pow(at2+m*m,2)));
    double brackets = m*m+(pplus*pplus+qplus*qplus)*at2/pow(pplus+qplus,2);
    return prefactor*brackets;
  }
  
  double qqg(const double p, const double phip, const double q, const double phiq, 
             const double k1, const double phik1, const double k2, const double phik2, 
             const double k, const double phik, const double yp, const double yq, const double m) {
    
    double pplus = sqrt(m*m+p*p)*exp(yp)/sqrt(2.);
    double pminus = sqrt(m*m+p*p)*exp(-yp)/sqrt(2.);
    double qplus = sqrt(m*m+q*q)*exp(yq)/sqrt(2.);
    double qminus = sqrt(m*m+q*q)*exp(-yq)/sqrt(2.);

    // a_T = q_T-k_T (2D vectors)
    double at2 = q*q + k*k - 2*q*k*cos(phiq-phik);
    // b_T = q_T-k_T-k1_T (2D vectors)
    double bt2 = at2 + k1*k1 -2.*q*k1*cos(phiq-phik1) + 2.*k*k1*cos(phik-phik1);
    double adotb = at2 + k*k1*cos(phik-phik1) - q*k1*cos(phiq-phik1);
    double adotq = q*q - k*q*cos(phik-phiq);
    double bdotq = adotq - k1*q*cos(phik1-phiq);
    double adotp = q*p*cos(phiq-phip) - k*p*cos(phik-phip);
    double bdotp = adotp - k1*p*cos(phik1-phip);

    //Lipatov vertex
    double Cplus = pplus + qplus - k1*k1/(pminus+qminus);
    double Cminus = k2*k2/(pplus+qplus) - pminus - qminus;
    double adotC = q*k2*cos(phiq-phik2) -q*k1*cos(phiq-phik1)
      -k*k2*cos(phik-phik2) +k*k1*cos(phik-phik1);
    double bdotC = adotC - k1*k2*cos(phik1-phik2) + k1*k1;
    double Cp = (Cplus*pminus + Cminus*pplus) - k2*p*cos(phik2-phip) + k1*p*cos(phik1-phip);
    double Cq = (Cplus*qminus + Cminus*qplus) - k2*q*cos(phik2-phiq) + k1*q*cos(phik1-phiq);
    double pq = (pplus*qminus + pminus*qplus) - p*q*cos(phip-phiq);
    
    return 4./((m*m+pq)*(pplus*(m*m+at2)+qplus*(m*m+bt2)))*
      ((m*m+adotb)*(qplus*Cp+pplus*Cq-Cplus*(m*m+pq))+Cplus*((m*m+bdotq)*(m*m-adotp)-(m*m+adotq)*(m*m-bdotp))+pplus*(adotC*(m*m+bdotq)-bdotC*(m*m+adotq))+qplus*(adotC*(m*m-bdotp)-bdotC*(m*m-adotp)));
    
  }

  double gg(const double p, const double phip, const double q, const double phiq, 
            const double k1, const double phik1, const double k2, const double phik2, 
            const double k, const double phik, const double yp, const double yq, const double m){

    double pplus = sqrt(m*m+p*p)*exp(yp)/sqrt(2.);
    double pminus = sqrt(m*m+p*p)*exp(-yp)/sqrt(2.);
    double qplus = sqrt(m*m+q*q)*exp(yq)/sqrt(2.);
    double qminus = sqrt(m*m+q*q)*exp(-yq)/sqrt(2.);

    double Cplus = pplus + qplus - k1*k1/(pminus+qminus);
    double Cminus = k2*k2/(pplus+qplus) - pminus - qminus;

    double Cp = (Cplus*pminus + Cminus*pplus) - k2*p*cos(phik2-phip) + k1*p*cos(phik1-phip);
    double Cq = (Cplus*qminus + Cminus*qplus) - k2*q*cos(phik2-phiq) + k1*q*cos(phik1-phiq);
    double pq = (pplus*qminus + pminus*qplus) - p*q*cos(phip-phiq);
    double C2 = 2.*Cplus*Cminus - ( k2*k2+k1*k1-2.*k1*k2*cos(phik2-phik1) );
    
    return (2.*Cp*Cq - (m*m+pq)*C2)/pow((m*m+pq),2.);
    
  }
  

}
