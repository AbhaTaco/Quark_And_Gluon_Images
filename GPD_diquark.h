#ifndef GPD_DIQUARK_H
#define GPD_DIQUARK_H


const int     NGPD     = 2     ;


// parameters taken from the flavor paper

/// u

const double  apu      = 2000. ;   // ERBL

const double  Lambdau  = 1.018 ;
const double  Mxu      = 0.604 ;
const double  mu       = 0.420 ;

/// d

const double  apd      = 1000. ;   // ERBL

const double  Lambdad  = 0.860 ;
const double  Mxd      = 0.913 ;
const double  md       = 0.275 ;

/// g

const double  apg      = 1000. ;   // ERBL

const double  Lambdag  = 0.979 ;
const double  Mxg      = 0.726 ;



//

const double  alphaHu  = 0.210 ;
const double  alphapHu = 2.448 ;
const double  pHu      = 0.620 ;
const double  NormHu   = 2.043 ;

const double  ealphapHu = 0.0885 ;
const double  epHu      = 0.0725 ;
const double  eNormHu   = 0.0000 ;

const double  alphaEu  = 0.210 ;
const double  alphapEu = 2.835 ;
const double  pEu      = 0.969 ;
const double  NormEu   = 1.803 ;

const double  ealphapEu = 0.1460 ;
const double  epEu      = 0.3355;
const double  eNormEu   = 0.0000 ;

////////////

const double  alphaHd  = 0.0317 ;
const double  alphapHd = 2.209 ;
const double  pHd      = 0.658 ;
const double  NormHd   = 1.570 ; 

const double  ealphapHd = 0.1560 ;
const double  epHd      = 0.2570 ;
const double  eNormHd   = 0.0000 ;

const double  alphaEd  = 0.0317 ;
const double  alphapEd = 1.281 ;
const double  pEd      = 0.726 ;
const double  NormEd   = -2.780 ;

const double  ealphapEd = 3.176 ;
const double  epEd      = 1.543 ;
const double  eNormEd   = 0.000 ;

////////////

const double  alphaHg  = -0.622 ;
const double  alphapHg = 2.000 ;
const double  pHg      = 2.000 ;
const double  NormHg   = 1.467 ;

const double  ealphapHg = 0.10 ;
const double  epHg      = 0.05 ;
const double  eNormHg   = 0.228 ;

const double  alphaEg  = -0.622 ;
const double  alphapEg = 0.000 ;
const double  pEg      = 0.000 ;
const double  NormEg   = 0.034 ;

const double  ealphapEg = 1.212 ;
const double  epEg      = 1.197 ;
const double  eNormEg   = 0.05 ;


/* * */

double FDGLAPkTphi(int GPD, int FLAV, double kx, double ky, double xs, double xi, double tt) ;

double FDGLAPkT(int GPD, int FLAV, double kT, double xs, double xi, double tt) ;

double FDGLAP(int GPD, int FLAV, double xs, double xi, double tt) ;

double eFDGLAPkT(int GPD, int FLAV, double kT, double xs, double xi, double tt) ;

double eFDGLAP(int GPD, int FLAV, double xs, double xi, double tt) ;

double FERBL(int GPD, int FLAV, double xs, double xi, double tt, int j) ;

double F(int GPD, int FLAV, double xs, double xi, double tt) ;

double eF(int GPD, int FLAV, double xs, double xi, double tt) ;

double Fmom(int GPD, int FLAV, int n, double xi, double tt) ;

 
/***/

#endif
