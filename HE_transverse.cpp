#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#include "util.h"
#include "FFTDeltaT.h"
#include "GPD_diquark.h"


int main(int argc, char** argv){
  
  double kx, ky, bx, x, t, dx, xmin ;

  double total, Hfft, Efft, shift ;
  
  int Nx = 10 ;

  int FLAV = atof(argv[1]) ;

  xmin = 0.001 ;
  dx   = 0.09 ;
  x    = xmin ;
    
  double* A = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  double* B = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  double* C = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  double* D = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  double* E = (double*)malloc(sizeof(double) * Nfft * Nfft) ;
  double* F = (double*)malloc(sizeof(double) * Nfft * Nfft) ;

  for(int k = 0 ; k < Nx ; k++){
    
    for(int i = 0 ; i < Nfft ; i++){
      kx = (double)(i - Nfft / 2) / sqrt((double) Nfft) ;
    
      for(int j = 0 ; j < Nfft ; j++){
	ky = (double)(j - Nfft / 2) / sqrt((double) Nfft) ;
	
	t  = kx * kx + ky * ky ;
	
	A[i + Nfft * j] = FDGLAP(0, FLAV, x, 0., t) ;
	B[i + Nfft * j] = kx / M_ * FDGLAP(1, 0, x, 0., t) ;

	C[i + Nfft * j] = FDGLAP(0, FLAV, x, 0., t) ;
	D[i + Nfft * j] = 0. ;

	E[i + Nfft * j] = FDGLAP(1, FLAV, x, 0., t) ;
	F[i + Nfft * j] = 0. ;
      }
    } 

    DTtobTfftc(A, B) ;

    DTtobTfftc(C, D) ;

    DTtobTfftc(E, F) ;


    for(int i = 0 ; i < Nfft ; i++){
      bx  = 2. * pi * (double)(i - Nfft / 2) / sqrt((double) Nfft) ;
      
      total = A[i + Nfft * Nfft / 2] ;

      Hfft  = C[i + Nfft * Nfft / 2] ; 

      Efft  = E[i + Nfft * Nfft / 2] ;
      
      shift = total - Hfft ;
      
      printf("%e,%e,%e,%e,%e,%e\n", bx, x, total, shift, Hfft, Efft) ;
    }

    printf("\n\n\n") ;

    x += dx ;
  }

  free((void*)A) ;
  free((void*)B) ;
  free((void*)C) ;
  free((void*)D) ;
  free((void*)E) ;
  free((void*)F) ;
    
  return 0 ;
 
}
