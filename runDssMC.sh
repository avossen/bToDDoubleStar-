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
#path add_module main gibbetnit
path add_condition main >:0:analysis
path add_condition main <=:0:KILL

module put_parameter bToDDoubleStar rfname\/people/voss771/bToDDoubleStar/mcEx55.root
module put_parameter fix_mdst Make_pi0_option\2
module put_parameter fix_mdst Make_pi0_lower_limit\0.12
module put_parameter fix_mdst Make_pi0_upper_limit\0.15
module put_parameter fix_mdst Correct_ecl_option\1
module put_parameter fix_mdst Correct_ecl_pv\1
module put_parameter fix_mdst Make_pi0_pv\1


initialize
nprocess set 0
histogram define ${HBKFILE}
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-1_fullrecon.mdst 0
file 1
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-2_fullrecon.mdst 0
file 2
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-3_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-4_fullrecon.mdst 0
file 3
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-5_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-6_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-7_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-8_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-9_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-10_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-11_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-12_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-13_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-14_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-15_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-16_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-17_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-18_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-19_fullrecon.mdst 0
file 10
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-20_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-21_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-22_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-23_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-24_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-25_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-26_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-27_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-28_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-29_fullrecon.mdst 0
process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/DssMC_EKPFULL/evtgen_exp_55_Dss-30_fullrecon.mdst 0


output close
terminate
EOF
