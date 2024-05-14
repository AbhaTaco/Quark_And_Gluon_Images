# GluonImaging

This project uses C++ code and Bash scripts. 

The code takes discrete Fourier transforms (DFT) of experimental data supplemented by theoretical models to create spatial density distributions of elementary particles inside the proton. The outputs are datafiles that store the density distributions as well as plots of the distribution vs spatial length, which is at the Femtometer scale. The Discrete Fourier transforms are performed using the FFTW library.    

GPDFourierPlot.sh: produces the spatial density distribution, for a particular kind of particle in a particular spin Orientation of the proton. It involves the following C++ files, 
- GPDFourier.cpp: performs the DFT on the theoretical model and experimental data.
- FFTDeltaT.cpp: contains functions that perform the DFT.
- GPD_diquark.cpp: theoretical model and experimental data.
- util.cpp: constants and other miscellaneous functions necessary to execute the code.



 
