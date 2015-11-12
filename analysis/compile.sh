#!/bin/bash
CFLAGS=" -Wall -ggdb `root-config --cflags --libs`  -lMinuit"
echo cflags: $CFLAGS
c++ $CFLAGS doAnalysis.cc -o doAnalysis

