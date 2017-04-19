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
#needed for index files
#path add module ekpturbo
path add_condition main >:0:analysis
path add_condition main <=:0:KILL

module put_parameter bToDDoubleStar rfname\mcEx55.root
#module put_parameter fix_mdst Make_pi0_option\2
#module put_parameter fix_mdst Make_pi0_lower_limit\0.12
#module put_parameter fix_mdst Make_pi0_upper_limit\0.15
#module put_parameter fix_mdst Correct_ecl_option\1
initialize
nprocess set 0
histogram define ${HBKFILE}
process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000003r000040-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst


process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000041r000066-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000155r000168-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000191r000206-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0

process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000347r000371-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000585r000599-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0

process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000870r000881-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0

process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000913r000935-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0

process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000969r000976-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r001096r001136-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r001191r001209-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r001346r001355-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r001373r001416-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r001417r001431-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r001709r001726-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000422r000464-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000207r000227-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000372r000421-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s00/xlnu-e000055r000466r000491-s00-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
 
 
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000003r000040-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000041r000066-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000068r000123-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000124r000137-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000138r000154-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000155r000168-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000170r000189-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000191r000206-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000207r000227-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000228r000261-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000262r000313-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000315r000346-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000347r000371-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000372r000421-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000422r000464-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000466r000491-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000492r000508-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000509r000525-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000527r000564-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000565r000584-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000585r000599-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000600r000612-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000618r000642-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000644r000656-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000657r000669-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000670r000682-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000684r000695-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /ghi/fs01/belle/bdata2/users/robin/xlnu/mc/on_resonance/e000055/evtgen-mixed/s01/xlnu-e000055r000697r000720-s01-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
#

 
# 
# 
# process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/XLNU/on_resonance/e000021/evtgen-mixed/s12/xlnu-e000021r000164r000324-s12-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/XLNU/on_resonance/e000021/evtgen-mixed/s12/xlnu-e000021r000164r000163-s12-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# process_event /pic/projects/belle/ZMDST/MC/MC_SPECIAL/XLNU/on_resonance/e000021/evtgen-mixed/s12/xlnu-e000021r000164r000096-s12-evtgen-mixed-on_resonance-b20090127_0910.mdst 0
# 
output close
terminate
EOF
