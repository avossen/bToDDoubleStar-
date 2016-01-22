#!/bin/sh
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=.:`root-config --libdir`:$LD_LIBRARY_PATH
export BASF_Module_Dir=.:$basf_Module_DIR
#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH
echo "Module dir..."
echo    $BASF_MODULE_DIR
echo "panther tables..."
echo    $PANTHER_TABLE_DIR

MODULE=bToDDoubleStar
HBKFILE=/dev/null

basf <<EOF
path create main
path create analysis
path add_module main fix_mdst
path add_module main ${MODULE}
path add_module main gibbetnit
path add_condition main >:0:analysis
path add_condition main <=:0:KILL

module put_parameter bToDDoubleStar rfname\/people/voss771/bToDDoubleStar/mcEx55.root
module put_parameter fix_mdst Make_pi0_option\2
module put_parameter fix_mdst Make_pi0_lower_limit\0.12
module put_parameter fix_mdst Make_pi0_upper_limit\0.15
module put_parameter fix_mdst Correct_ecl_option\1
initialize
nprocess set 0
histogram define ${HBKFILE}

process_event /pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/on_resonance/mixed/evtgen-mixed-00-all-e000055r000065-b20090127_0910.mdst 0  

output close
terminate
EOF
