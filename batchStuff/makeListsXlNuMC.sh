#!/bin/bash
#example
# /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-charged/s01/
for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for spec in mixed charged uds charm
do
for stream in 0 1 2 3 4 
do
find /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e0000$ex/evtgen-$spec/s0$stream/  -iname '*.mdst' > listsXlNuMC/xlnuMC_ex$ex\_$res\_$spec\_$stream.list
find /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e0000$ex/evtgen-$spec/s1$stream/  -iname '*.mdst' >> listsXlNuMC/xlnuMC_ex$ex\_$res\_$spec\_$stream.list
done
done
done

#/pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/continuum/charm/