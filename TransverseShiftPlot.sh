#!/bin/bash


#
# bash TransverseShiftPlot.sh FLAV
#
# FLAV 0 u 1 d 2 g
#

# # #


filnm=HE_transverse

narg=$#


# * * * * * * * * * *


declare -a sFLAV_ar=(u d g)

FLAV=$1

sFLAV=${sFLAV_ar[$FLAV]}


# * * * * * * * * * *


home_dir=$HOME/Documents/Research/gluon_imaging

main_dir=${home_dir}/FourierTransforms/FourierOfGPDs

Fourier_dir=${home_dir}/FourierTransforms

kelley_dir=${home_dir}/kelley

export CPLUS_INCLUDE_PATH=${kelley_dir}:${Fourier_dir}:$main_dir:$home_dir:$CPLUS_INCLUDE_PATH


# * * * * * * * * * * 


main=${main_dir}/${filnm}.cpp

imagefile=${main_dir}/data_plots/${sFLAV}Transverse.eps

datafile=${main_dir}/data_plots/${sGPD}${sFLAV}Transverse.dat

g++ -w ${home_dir}/util.cpp ${kelley_dir}/FF_kelley.cpp ${home_dir}/GPD_diquark.cpp ${Fourier_dir}/FFTDeltaT.cpp $main -lfftw3 -o go

#./go $FLAV > $datafile

printf "$narg\n"

# * * * * * * * * * * * * * * * * * * * * 
#               GNUPLOT                 #
# * * * * * * * * * * * * * * * * * * * *  


#if [[ ${narg} < 10 ]] ; 

   
gnuplot <<EOF

set term postscript eps enhanced defaultplex \
    leveldefault color colortext \
    size 6.50in, 4.00in "Helvetica" 15 fontscale 1.5

set output '$imagefile'

set datafile separator ','

set xlabel 'b_T (GeV^{-1})'

#set title "${sGPD}"."${sFLAV}".'(x, 0, b_T)'

set xrange [-20 : 20]
#set yrange [0 : ]

set label "${sFLAV}" at screen 0.3,0.7

plot '$datafile' index 0 u 1:5 t "x 0.001"  with lines lc 0 lw 3,\
'' index 0 u 1:3 notitle with lines dt 2 lc 0 lw 3,\
'' index 3 u 1:5 t "x 0.27" with lines lc 1 lw 3,\
'' index 3 u 1:3 notitle with lines dt 2 lc 1 lw 3,\
'' index 6 u 1:5 t "x 0.54" with lines lc 3 lw 3,\
'' index 6 u 1:3 notitle with lines dt 2 lc 3 lw 3
#'' index 8 u 1:5 t "x 0.72" with lines lc 4 lw 3,\
#'' index 8 u 1:3 notitle with lines dt 2 lc 4 lw 3


EOF


rm *~

#fi

# # #
