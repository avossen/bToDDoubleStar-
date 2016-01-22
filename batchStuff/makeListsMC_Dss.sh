#!/bin/bash

for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
find /pic/projects/Belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/ -iname "*exp_$ex*.mdst" > lists/mc$ex.list
done


#/pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/continuum/charm/