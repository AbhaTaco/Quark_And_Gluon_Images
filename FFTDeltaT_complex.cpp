#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#include "util.h"
#include "FFTDeltaT_complex.h"




void DTtobTfftc(double* Re, double* Im){

  double norm ;

  norm = 2. * pi * 2. * pi * (double) Nfft ;
  
  fftw_plan p ;

  fftw_complex*  in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft * Nfft) ;

  for(int i = 0; i < Nfft; i++)
    for(int j = 0; j < Nfft; j++){
      in[i + Nfft * j][0] = Re[i + Nfft * j] ;
      in[i + Nfft * j][1] = Im[i + Nfft * j] ;
    }
  
  p = fftw_plan_dft_2d(Nfft, Nfft, in, in, FFTW_FORWARD, FFTW_ESTIMATE) ;  

  fftw_execute(p) ;   

  for(int i = 0; i < Nfft; i++)
    for(int j = 0; j < Nfft; j++){
      Re[i + Nfft * j] = in[(i + Nfft/2) % Nfft + Nfft * ((j + Nfft/2) % Nfft)][0] * pow(-1, i + j) / norm ;
      Im[i + Nfft * j] = in[(i + Nfft/2) % Nfft + Nfft * ((j + Nfft/2) % Nfft)][1] * pow(-1, i + j) / norm ;
    }

  fftw_destroy_plan(p) ;
  fftw_free(in); 

  
}



/* * * * * * * * * * * * * * * * * * */



void bTtoDTfftc(double* Re, double* Im){

  double norm ;

  norm = (double) Nfft / (4. * pi * pi) ;
  
  fftw_plan p ;

  fftw_complex*  in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft * Nfft) ;

  for(int i = 0; i < Nfft; i++)
    for(int j = 0; j < Nfft; j++){
      in[i + Nfft * j][0] = Re[i + Nfft * j] ;
      in[i + Nfft * j][1] = Im[i + Nfft * j] ;
    }
  
  p = fftw_plan_dft_2d(Nfft, Nfft, in, in, FFTW_BACKWARD, FFTW_ESTIMATE) ;  

  fftw_execute(p) ;   

  for(int i = 0; i < Nfft; i++)
    for(int j = 0; j < Nfft; j++){
      Re[i + Nfft * j] = in[(i + Nfft/2) % Nfft + Nfft * ((j + Nfft/2) % Nfft)][0] * pow(-1, i + j) / norm ;
      Im[i + Nfft * j] = in[(i + Nfft/2) % Nfft + Nfft * ((j + Nfft/2) % Nfft)][1] * pow(-1, i + j) / norm ;
    }

  fftw_destroy_plan(p) ;
  fftw_free(in); 

}



/* * * * * * * * * * * * * * * * * * */
