
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#include "util.h"
#include "FFTDeltaT.h"
#include "GPD_diquark.h"


int main(int argc, char** argv){
  
    
  double kx, ky, bx, cal, x, t, dx, xmin, eA, ecal ;

  int Nx = 8 ;

  int nar = argc - 1 ;

  int GPD = atof(argv[1]) ;
  int FLAV = atof(argv[2]) ;
  
  
  xmin = 0.1 ;
  dx   = 0.1 ;
  x    = xmin ;
    
  double* A = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  double* Au = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  double* Al = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  
  for(int k = 0 ; k < Nx ; k++){
    
    for(int i = 0 ; i < Nfft ; i++){
      kx = (double)i / 2. / sqrt((double) Nfft) ;
    
      for(int j = 0 ; j < Nfft ; j++){
	ky = (double)j / 2. / sqrt((double) Nfft) ;

	t = kx * kx + ky * ky ;

	eA = eFDGLAP(GPD, FLAV, x, 0., t) ;
	
	A[i + Nfft * j] = FDGLAP(GPD, FLAV, x, 0., t) ;

	Au[i + Nfft * j] = A[i + Nfft * j] + eA ;

	Al[i + Nfft * j] = A[i + Nfft * j] - eA ;        

      }
    }

    
    DTtobTfft(A) ;

    DTtobTfft(Au) ;

    DTtobTfft(Al) ;

    
    for(int i = 0 ; i < Nfft ; i++){
      bx = 2. * pi * (double)i / sqrt((double) Nfft) ;
      
      cal = A[i] ;

      ecal = (fabs(A[i] - Au[i]) + fabs(A[i] - Al[i])) / 2. ;
	
      printf("%e,%e,%e,%e,%e\n", bx, x, cal, cal + ecal, cal - ecal) ;
    }

    //    printf("\n\n\n") ;

    x += dx ;
  }

  free((void*) A) ;

  free((void*) Au) ;

  free((void*) Al) ;
  
  return 0 ;
 
}
