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
process_event Simulation/mdst/evtgen_exp_07_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_09_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_11_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_13_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_15_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_17_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_19a_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_19b_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_21_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_23_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_25_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_27_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_31_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_33_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_35_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_37a_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_37b_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_39_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_41a_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_41b_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_43_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_45a_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_45b_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_47_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_49_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_51_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_55_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_55_bToDPiPi-1.mdst
process_event Simulation/mdst/evtgen_exp_61_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_63_bToDPiPi-0.mdst
process_event Simulation/mdst/evtgen_exp_65_bToDPiPi-0.mdst


output close
terminate
EOF
