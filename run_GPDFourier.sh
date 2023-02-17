
#!/bin/bash

#
# bash run_GPDFourier.sh
#

# # #



for(( i=0 ; i < 2 ; i++))
    do
	for(( j=0 ; j < 3 ; j++))
	do

	    bash GPDFourierPlot.sh $i $j
	    
	done
done


#GPD=0 FLAV=0
#bash GPDFourierPlot.sh $GPD $FLAV

# # #
