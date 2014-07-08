#!/bin/bash

DIR=`echo $1`
ENERGY=`echo $2`

echo ""
echo ">> `/bin/date` Submitting Condor job(s)..."
#==========================================================
#==========================================================
#==========================================================
echo ">> `/bin/date` PythiaMonoJet: parameters in $1"

mkdir ${DIR} -p
cat ./PythiaMonoJet_cfi_tmp.py \
| sed -e s/ENERGYIN/${ENERGY}/ \
> ./${DIR}/PythiaMonoJet_cfi_${ENERGY}.py

cat ./FScheck.jdl \
| sed -e s/PREFIX_NAME/${DIR}/ \
| sed -e s~RUN_DIR~/home/pedrok/CMSSW_4_2_8/src/ForwardCaloUpgrade/MCtest/test/${DIR}~ \
| sed -e s/INPUT_NAME/PythiaMonoJet_cfi_${ENERGY}.py/ \
> FS_MonoJet_${ENERGY}.jdl

condor_submit FS_MonoJet_${ENERGY}.jdl
