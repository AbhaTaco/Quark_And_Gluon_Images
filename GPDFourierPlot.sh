#!/bin/bash


#
# bash GPDFourierPlot.sh GPD FLAV
#
# GPD 0 H 1 E
# FLAV 0 u 1 d 2 g
#

# # #


filnm=GPDFourier

narg=$#


# * * * * * * * * * *


declare -a sGPD_ar=(H E)
declare -a sFLAV_ar=(u d g)

GPD=$1

FLAV=$2

sGPD=${sGPD_ar[$GPD]}
sFLAV=${sFLAV_ar[$FLAV]}


# * * * * * * * * * *


home_dir=$HOME/Documents/Research/gluon_imaging

main_dir=${home_dir}/FourierTransforms/FourierOfGPDs

Fourier_dir=${home_dir}/FourierTransforms

kelley_dir=${home_dir}/kelley

export CPLUS_INCLUDE_PATH=${kelley_dir}:${Fourier_dir}:$main_dir:$home_dir:$CPLUS_INCLUDE_PATH


# * * * * * * * * * * 


main=${main_dir}/${filnm}.cpp

imagefile=${main_dir}/data_plots/${sGPD}${sFLAV}Fourier3D.pdf

datafile=${main_dir}/data_plots/${sGPD}${sFLAV}Fourier3D.dat

g++ -w ${home_dir}/util.cpp ${kelley_dir}/FF_kelley.cpp ${home_dir}/GPD_diquark.cpp ${Fourier_dir}/FFTDeltaT.cpp $main -lfftw3 -o go

#./go $GPD $FLAV > $datafile

printf "$narg\n"

# * * * * * * * * * * * * * * * * * * * * 
#               GNUPLOT                 #
# * * * * * * * * * * * * * * * * * * * *  


#if [[ ${narg} > 10 ]] ; then 

   
gnuplot <<EOF


set term pdf font "Helvetica, 15"  \
    size  6.50in, 5.00in 

set output '$imagefile'

set datafile separator ','

set style fill transparent solid 0.5 border
set style data filledcurves 

set xrange [0:5]

set yrange [0.3:0.6]

#set zrange [0:]

set xyplane at 0
# rotate by -12
set xlabel "b_T (GeV^{-1})" offset -1,-1.5

set border 895 #319 #63

set grid vertical 
set grid ztics

set ylabel "x" offset 1,0

set xtics nomirror offset 0,-0.5

set ytics 0.1 scale 1 

#set ztics 0.05

set view 63,30

splot '$datafile' u 1:2:3 t "H_g" w l lw 3 lc 'blue',\
'' u 1:2:4:5 notitle 

EOF

#fi

# # #
