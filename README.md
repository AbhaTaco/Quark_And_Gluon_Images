# GluonImaging

This project uses C++ code and Bash scripts. 

The code takes discrete Fourier transforms (DFT) of experimental data supplemented by theoretical models to create spatial density distributions of elementary particles inside the proton. The outputs are datafiles that store the density distributions as well as plots of the distribution vs spatial length, which is at the Femtometer scale. The Discrete Fourier transforms are performed using the FFTW library.    

**GPDFourierPlot.sh**: produces the spatial density distribution, for a particular kind of particle (u,d or g depicted by 0,1 or 2 repectively) in a particular spin orientation of the proton, H and E depicted by 0 or 1. It involves the following C++ files, 
- **GPDFourier.cpp**: performs the DFT on the theoretical model and experimental data.
- **FFTDeltaT.cpp**: contains functions that perform the DFT.
- **GPD_diquark.cpp**: theoretical model and experimental data.
- **util.cpp**: constants and other miscellaneous functions necessary to execute the code.

**TransverseShiftPlots.sh**: produces the spatial density distribution, but this time for a transversely polarized proton 
- **HE_transverse.cpp**: performs the DFT on the theoretical model and experimental data for a transversely polarized proton.  
The rest of the files are the same as used in GPDFourierPlot.sh. 


 
