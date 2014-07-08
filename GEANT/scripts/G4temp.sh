#!/bin/bash

FILE=`echo $1`
INPUT=`echo $2`

echo ""
echo ">> `/bin/date` Submitting Condor job(s)..."
#==========================================================
#==========================================================
#==========================================================
echo ">> `/bin/date` fcalor: parameters in $2"

cat ./G4check.jdl \
| sed -e s/PREFIX_NAME/${FILE}/ \
| sed -e s~RUN_DIR~/home/pedrok/CMSSW_4_2_8/src/ForwardCaloUpgrade/GEANT-00_03_00~ \
| sed -e s/INPUT_NAME/${INPUT}/ \
> G4_${FILE}.jdl

condor_submit G4_${FILE}.jdl
