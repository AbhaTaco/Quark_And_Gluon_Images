#ifndef UTIL_H
#define UTIL_H


const double  M_     = 0.938 ;
const double  pi     = 4. * atan(1.) ;


/***/

void gauss(int npts, int job, double a, double b, double* x, double* w) ;

double calM (double x, double Mx, double Lambda) ;

double regf (double xs, double xi, double tt, double alpha, double alphap, double p) ;

double DTf (double xi, double tt) ;

double Den0 (double kT, double xs, double xi, double tt, double Mx, double Lambda) ;

double Den1 (double kT, double xs, double xi, double tt, double Mx, double Lambda) ;


#endif
