#ifndef FFTDELTAT_H
#define FFTDELTAT_H



const int Nfft = 500 ;



/***/

void DTtobTfft(double* A) ;

void bTtoDTfft(double* A) ;

void DTtobTfftc(double* A, double* B) ;

void bTtoDTfftc(double* A, double* B) ;


/***/


#endif
