#!/bin/bash

FILE=`echo $1`
INPUT=`echo $2`
ENERGY=`echo $3`
PART=`echo $4`

echo ""
echo ">> `/bin/date` Submitting Condor job(s)..."
#==========================================================
#==========================================================
#==========================================================
echo ">> `/bin/date` fcalor: parameters in $2"

cat ./${INPUT} \
| sed -e s/ENERGYIN/${ENERGY}/ \
| sed -e s/PNUMBER/${PART}/ \
> ${FILE}_${ENERGY}_${PART}.in

cat ./G4check.jdl \
| sed -e s/PREFIX_NAME/${FILE}/ \
| sed -e s~RUN_DIR~/home/pedrok/CMSSW_4_2_8/src/ForwardCaloUpgrade/GEANT-00_03_00~ \
| sed -e s/INPUT_NAME/${FILE}_${ENERGY}_${PART}.in/ \
> G4_${FILE}_${ENERGY}_${PART}.jdl

condor_submit G4_${FILE}_${ENERGY}_${PART}.jdl
