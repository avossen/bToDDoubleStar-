#!/bin/bash

#use like so: ./doAnalysis.sh ../SumMC_mixed.root ../SumMC_charged.root ../SumMC_uds.root ../SumMC_charm.root  2
root.exe -b -q ./doAnalysis.C\(\"$1\",\"$2\",\"$3\",\"$4\",$5\)