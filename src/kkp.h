//
//  kkp.h
//  
//
//  Created by Bjoern Schenke on 4/6/20.
//

#ifndef kkp_h
#define kkp_h

#include <stdio.h>
#include <iostream>             
#include <math.h>       

#endif /* kkp_h */

//#include "particle.h"

enum FLAVOR_TYPE {                                                                                                            
  gluon,                  ///< 0: g (gluon)                                                                                   
  up,                     ///< 1: u (up)                                                                                      
  anti_up,                ///< 2: ub (anti-up)                                                                                
  down,                   ///< 3: d (down)                                                                                    
  anti_down,              ///< 4: db (anti-down)                                                                              
  strange,                ///< 5: s (strange)                                                                                 
  anti_strange,           ///< 6: sb (anti-strange)                                                                           
  charm,                  ///< 7: c (charm)                                                                                   
  anti_charm,             ///< 8: cb (anti-charm)                                                                             
  bottom,                 ///< 9: b (bottom)                                                                                  
  anti_bottom,            ///< 10: bb (anti-bottom)                                                                           
  photon          =   13, ///< 13: p (photon)                                                                                 
  electron        = 1000, ///< 1000: electron                                                                                 
  positron        = 1001, ///< 1001: positron                                                                                 
  // generalized flavor type which should not assigned to an individual                                                       
  // particle, but is used for analysis purposes                                                                              
  quark           =   77, ///< 77: quark = light quark or anti-light-quark                                                    
  light_quark     =   88, ///< 88: light quark (generic)                                                                      
  anti_light_quark=   99, ///< 99: light anti-quark (generic)                                                                 
  allFlavors      =  111, ///< 111: any flavor (generic)                                                                      
  electron_gen    =  335       ///< 335:  all electrons/positrons (generic)                                                  
};                                                  


namespace kkp
    {
      double KKPFragmentation(int ih, int iset,   double x, double qs, FLAVOR_TYPE species);
    }   
