#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`


for stream in s00 s01 s02 s03 s04
do
for res in on_resonance continuum
do
for spec in mixed charged uds charm
do
for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do

myDir=subXlNuMC_ex$ex\_$res\_$spec\_$stream
mkdir $myDir
echo " dir: $myDir " ;
mkdir /pic/projects/belle/voss771/bToD/$myDir
mkdir /pic/projects/belle/voss771/bToD/$myDir


#for d in `cat mc$1_onRes.list`
for d in `cat listsFullMC/mc$ex\_$res\_$spec_\$stream.list`
do
let "counter+=1"
#if [ "$counter" -le "$5" ]
#then
#if [ "$counter" -ge "$4" ]
#then
let "subCounter+=1"
targetShFile=$myDir/job_$ex\_$res\_$spec\_$stream\_$subCounter.sh
#cp batchHead.sh $targetShFile
cp batchHead1.sh $targetShFile
echo "#SBATCH -o /pic/projects/belle/voss771/bToD/$myDir/jobId_$subCounter.out" >> $targetShFile
echo "#SBATCH -e /pic/projects/belle/voss771/bToD/$myDir/jobId_$subCounter.err" >> $targetShFile
echo "#SBATCH -J bToD_$subCounter"  >> $targetShFile 
cat batchHead2.sh >> $targetShFile 
echo "module put_parameter bToDDoubleStar rfname\\/pic/projects/belle/voss771/bToD/$myDir/job_$subCounter.root" >> $targetShFile
cat batchMiddle.sh >> $targetShFile
echo "process_event $d 0" >> $targetShFile
cat batchEnd.sh >> $targetShFile
#fi
#fi
done
done
done
done
done




