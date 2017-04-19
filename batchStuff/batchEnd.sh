
output close

terminate

EOF

dateString=`date +%d%b%Y`

# Grab the exit code of BASF
BASFRET=$?
echo VOSSEN_BASF_FINISH `date`


# Mark the file as bad if BASF returned something other than 0
if [ $BASFRET -ne 0 ]; then
  for p in /group/belle/users/vossen/bToDOut/*${SLURM_JOBID}*; do
    mv $p $p.badret
  done
fi

