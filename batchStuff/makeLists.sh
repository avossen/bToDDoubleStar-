#!/bin/bash

for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for kind in charged mixed
do
for period in p1 p2 
do
echo finding '*exp_'$ex'*'$kind'*'$period'.mdst'
find /group/belle/users_old/jansh/DssMC2/mdst/ -iname '*exp_'$ex'*'$kind'*'$period'*.mdst' > newHuschleMC$ex\_$kind\_$period.list
done
done
done