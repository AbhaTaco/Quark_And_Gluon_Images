#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "GPD_diquark.h"
#include "util.h"
#include "FF_kelley.h"

using namespace std ;

// t = - Delta * Delta in this code


const double m      [2] = {mu, md} ;
const double Mx     [3] = {Mxu, Mxd, Mxg} ;
const double Lambda [3] = {Lambdau, Lambdad, Lambdag} ;
const double ap     [3] = {apu, apd, apg} ;

const double Norm   [6] = {NormHu, NormEu, NormHd, NormEd, NormHg, NormEg} ;
const double alpha  [6] = {alphaHu, alphaEu, alphaHd, alphaEd, alphaHg, alphaEg} ;
const double alphap [6] = {alphapHu, alphapEu, alphapHd, alphapEd, alphapHg, alphapEg} ;
const double p      [6] = {pHu, pEu, pHd, pEd, pHg, pEg} ;

const double eNorm   [6] = {eNormHu, eNormEu, eNormHd, eNormEd, eNormHg, eNormEg} ;
const double ealphap [6] = {ealphapHu, ealphapEu, ealphapHd, ealphapEd, ealphapHg, ealphapEg} ;
const double ep      [6] = {epHu, epEu, epHd, epEd, epHg, epEg} ;


/* * * * */


double FDGLAPkTphi(int GPD, int FLAV, double kx, double ky, double xs, double xi, double tt){
  
  // GPD  -> 0 H, 1 E
  // FLAV -> 0 u, 1 d
  
  double a, Den, Denp, kT, x, zeta, DT, reg, xp ;

  int i = GPD + FLAV * NGPD ;
  int j = FLAV ;
  
  a = Den = Denp = 0. ;

  
  /*   kinematics   */
    
  kT = sqrt(kx*kx + ky*ky) ;

  zeta = 2. * xi / (1. + xi) ;
  x    = xs * (1. - zeta / 2.) + zeta / 2. ;
    
  xp   = (x - zeta) / (1. - zeta) ;
  
  DT   = sqrt( tt * (1. - zeta) - (M_ * M_ * zeta * zeta) ) ;

  /* * * * */

  
  reg  = pow( x, -alpha[i] + alphap[i] * pow(1. - x, p[i]) * tt) ;

  Den  = calM(x, Mx[j], Lambda[j]) - kT * kT / (1. - x) ;

  Denp = calM(xp, Mx[j], Lambda[j]) - kT * kT / (1. - xp) - (1. - xp) * DT * DT  + 2. * kx * DT ;

  if (GPD == 0){
    
    if (zeta > 0.)
      a = reg / (1. - x) * Norm[i] * (1. - zeta / 2.) * ((m[j] + x * M_) * (m[j] + xp * M_) + kT * kT - (1. - xp) * kx * DT ) / (Den * Den * Denp * Denp) +  zeta * zeta * FDGLAPkTphi(1, FLAV, kx, ky, xs, xi, tt) / (4. * (1. - zeta)) ;
    else
      a = reg / (1. - x) * Norm[i] * ((m[j] + x * M_) * (m[j] + x * M_) + kT * kT - (1. - x) * kx * DT ) / (Den * Den * Denp * Denp) ;
  }
  else if (GPD == 1)
      a    = - 2. * reg * M_ * Norm[i] * (1. - zeta / 2.) * (1. - zeta) * (- (M_ + x * m[j]) * (1. - xp) + (x - xp) * M_ * kx / DT) / (Den * Den * Denp * Denp) ;

  return a ;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



double FDGLAPkT(int GPD, int FLAV, double kT, double xs, double xi, double tt){
  
  // GPD  -> 0 H, 1 E
  // FLAV -> 0 u, 1 d, 2 g
  
  double an, Den, Denp, x, zeta, DT, reg, xp, a, b, beta ;

  int i = GPD + FLAV * NGPD ;
  int j = FLAV ;
  
  an = Den = Denp = 0. ;

  beta = 10. ;
  
  /*   kinematics   */
    
  zeta = 2. * xi / (1. + xi) ;
  x    = xs * (1. - zeta / 2.) + zeta / 2. ;
    
  xp   = (x - zeta) / (1. - zeta) ;
  
  DT   = sqrt( tt * (1. - zeta) - (M_ * M_ * zeta * zeta) ) ;

  /* * * * */

  
  reg  = pow( x, -alpha[i] + alphap[i] * pow(1. - x, p[i]) * tt)  ; //+ beta * tt * zeta * zeta / (1. - zeta)) ;

  a    = calM(xp, Mx[j], Lambda[j]) - kT * kT / (1. - xp) - (1. - xp) * DT * DT ; //took away overall minus sign Jan 29, 2021
    
  b    = 2. * DT * kT ; // b has to be positive so what we are using here is actually (-(kpsq - Lambda^2))^2

  Den  = calM(x, Mx[j], Lambda[j]) - kT * kT / (1. - x) ;

  Denp = (fabs(a * a - b * b)) * sqrt(fabs(a * a - b * b)) ;

  
  ////////

  if (FLAV == 2){

    
    if (GPD == 0)

      an   = - 2. * pi * reg / (1. - x) * Norm[i] * (a * ( x * xp * ((1. - x) * M_  - Mx[j]) * ((1. - xp) * M_  - Mx[j]) + (1. / (1. - xp) - (1. - x)) * kT * kT) - b * ( 1. / (1. - x) + 1. - xp) * kT * DT) / (Den * Den * Denp) + zeta * zeta / (4. * (1. - zeta)) * FDGLAPkT (1, FLAV, kT, xs, xi, tt) ;
  
    else if (GPD == 1)
    
      an = - 2. * pi * reg * Norm[i] * 2. * M_ * (1. - zeta) / (1. - zeta / 2. ) * ( 2. * kT * kT * (x * ((1. - x) * M_ - Mx[j]) - xp * (1. - zeta) * ((1. - xp) * M_ - Mx[j])) - a * (1. - xp) * x * ((1. - x) * M_ - Mx[j])) / (Den * Den * Denp) ;

    
  }
  
  else {
      
      
    if (GPD == 0)

      an   = -2. * pi * reg / (1. - x) * Norm[i] * (1. - zeta / 2.) * (a * ((m[j] + x * M_) * (m[j] + xp * M_) + kT * kT) - b * (1. -xp) * kT * DT) / (Den * Den * Denp) + zeta * zeta / (4. * (1. - zeta)) * FDGLAPkT (1, FLAV, kT, xs, xi, tt) ;
  
    else if (GPD == 1)
    
      an = -2. * pi * reg / (1. - x) * Norm[i] * (1. - zeta / 2.) / (1. - zeta) * (-4. * M_ * M_ * kT * kT * (x - xp) + a * 2. * M_ * (m[j] + x * M_) * (1. - xp)) / (Den * Den * Denp) ;
  
  }

  
  return an ;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



double FDGLAP(int GPD, int FLAV, double xs, double xi, double tt){

  double a, kTmax, kTmin, pp ;

  int Ng = 51 ;
  
  kTmin = 0. ;
  kTmax = 10. ;
  
  a     = 0.  ; 
    
  double* kT = (double*) calloc(Ng, sizeof(double)) ; 
  double* w = (double*) calloc(Ng, sizeof(double)) ;

  gauss(Ng, 0, kTmin, kTmax, kT, w) ;

  
  for(int i = 0 ; i < Ng ; i++){
    
    pp = FDGLAPkT(GPD, FLAV, kT[i], xs, xi, tt) ;
    a += w[i] * kT[i] * pp ;
  }
 
  free((void*)kT) ;
  free((void*)w) ;

  if (xs < xi)
    a = 0. ;
  
  return a ;

}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


double eFDGLAPkT(int GPD, int FLAV, double kT, double xs, double xi, double tt){
  
  // GPD  -> 0 H, 1 E
  // FLAV -> 0 u, 1 d, 2 g
  
  double e1, e2, e3, e ;

  int i = GPD + FLAV * NGPD ;
  int j = FLAV ;

  double zeta = 2. * xi / (1. + xi) ;
  
  double x = xs * (1. - zeta / 2.) + zeta / 2. ;
  
  
  e1 = - pow((1. - x), p[i]) * tt * log(x) * ealphap[i] ; // the error from alpha

  e2 = - pow((1. - x), p[i]) * tt * alphap[i] * log(1. - x) * log(x) * ep[i] ; // the error from p

  e3 =  eNorm[i] / Norm[i] ;
    
  e  =  FDGLAPkT(GPD, FLAV, kT, xs, xi, tt) * sqrt (e1 * e1 + e2 * e2 + e3 * e3) ;

  
  return e ;
  
}


/* * * * * * * * * * * * * * * */


double eFDGLAP(int GPD, int FLAV, double xs, double xi, double tt){

  double a, kTmax, kTmin, pp ;

  int Ng = 51 ;
  
  kTmin = 0. ;
  kTmax = 10. ;
  
  a     = 0.  ; 
    
  double* kT = (double*) calloc(Ng, sizeof(double)) ; 
  double* w = (double*) calloc(Ng, sizeof(double)) ;

  gauss(Ng, 0, kTmin, kTmax, kT, w) ;

  
  for(int i = 0 ; i < Ng ; i++){
    
    pp = eFDGLAPkT(GPD, FLAV, kT[i], xs, xi, tt) ;
    a += w[i] * kT[i] * pp * w[i] * kT[i] * pp ;
  }

  a = sqrt(a) ;
  
  free((void*)kT) ;
  free((void*)w) ;

  if (xs < xi)
    a = 0. ;
  
  return a ;

}


/* * * * * * * * * * * * * * * */


double FERBL(int GPD, int FLAV, double xs, double xi, double tt, int j){

  double a, S, A, zeta, x, xmin, xmax, pp ;

  double am, a1p, a2p, Fm, Fp, az, c, d, an ;
  
  zeta = 2. * xi / (1. + xi) ;
  x    = xs * (1. - zeta / 2.) + zeta / 2. ;

  az   = FDGLAP(GPD, FLAV, xi, xi, tt) ;
  
  a1p  = ap [FLAV] ;
  a2p  = az * 2. / zeta - a1p * zeta * zeta / 4. ;

  c    = ( 2. * az + 0.5 * a1p * zeta * zeta * zeta ) / zeta ;
  d    = - az ;
 
  int Ng = 51 ;

  xmin  = xi    ;
  xmax  = 1.    ;
  A     = 0.    ; 
    
  double* xv  = (double*) calloc(Ng, sizeof(double)) ; 
  double* w   = (double*) calloc(Ng, sizeof(double)) ;

  gauss(Ng, 0, xmin, xmax, xv, w) ;

  
  for(int i = 0 ; i < Ng ; i++){
    
    pp  = FDGLAP(GPD, FLAV, xv[i], xi, tt) ;
    A  += w[i] * pp ;
  }
  
  S  = (1. - zeta / 2.) * (FF(GPD, FLAV, tt) - A) ;

  am = 6. * (zeta * az - 2. * S) / (zeta * zeta * zeta) ;

  Fm = am * x * (x - zeta) + az ;

  //  Fp = a1p * (x - zeta/2.) * (x - zeta/2.) * (x - zeta/2.) + a2p  * (x - zeta/2.) ;

  Fp = a1p * x * x * (x - 1.5 * zeta) + c * x + d ;
    
  a  = (Fm + Fp) / 2. ;

  // printf("%e,%e,%e,%e,%e,%e\n", xs, xi, tt, fabs(xs) * Fm, fabs(xs) * Fp, fabs(xs) * a) ;

  //printf("%e,%e,%e,%e,%e,%e\n", xs, xi, tt, xs * Fm, xs * Fp, xs * a) ;
  
  // printf("%e,%e,%e,%e,%e,%e\n", x, zeta, tt, Fm, Fp, a) ;
  //  printf("%e,%e,%e,%e,%e,%e\n", x, zeta, tt, fabs(x) * Fm, fabs(x) * Fp, fabs(x)*a) ;
  // printf("%e,%e,%e,%e,%e,%e,%e\n", xs, xi, tt, S, az, a2p, am) ;
  
  free((void*)xv) ;
  free((void*)w) ;

  if (j == 0)
    an = a ;
  else if (j == 1)
    an = Fm ;
  else if (j == 2)
    an = Fp ;
  
  if (xs > xi || xs == xi)
    an = 0. ;

  
  return an ;


}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


double F(int GPD, int FLAV, double xs, double xi, double tt){

  double a = 0. ;

  if( xs > 0. && xs < 1. )
    a = FDGLAP( GPD, FLAV, xs, xi, tt) ;
  else
    a = 0. ;
  

  return a ;
  
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


double eF(int GPD, int FLAV, double xs, double xi, double tt){

  double a = 0. ;

  if( xs > 0. && xs < 1. )
    a = eFDGLAP(GPD, FLAV, xs, xi, tt) ;
  else
    a = 0. ;


  return a ;
  
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



double Fmom(int GPD, int FLAV, int n, double xi, double tt){

  double a, xmax, xmin, pp ;

  int Ng = 51   ;

  xmin   = -xi  ;
  xmax   = 1.   ;
  a      = 0.   ; 
    
  double* x = (double*) calloc(Ng, sizeof(double)) ; 
  double* w  = (double*) calloc(Ng, sizeof(double)) ;

  gauss(Ng, 0, xmin, xmax, x, w) ;

  
  for(int i = 0 ; i < Ng ; i++){
    
    pp = pow(x[i], n) * F(GPD, FLAV, x[i], xi, tt) ;
    a += w[i] * pp ;
  }
 
  free((void*)x) ;
  free((void*)w) ;

  //a *= -1. ;
  
  return a ;


}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
