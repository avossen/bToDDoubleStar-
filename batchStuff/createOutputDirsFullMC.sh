#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`



for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
myDir=subXlNuMC_ex$ex
echo " dir: $myDir " ;
mkdir /pic/projects/belle/voss771/bToD/$myDir
mkdir /pic/projects/belle/voss771/bToDOut/$myDir

done








