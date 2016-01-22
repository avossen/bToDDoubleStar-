#!/bin/bash



for res in on_resonance continuum
do
for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
hadd -f SumXlNuData_Ex$ex\_$res.root /pic/projects/belle/voss771/bToD/subXlNuData_ex$ex\_$res/*.root
done
done
#/pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/continuum/charm/