#!/bin/bash

DIR=`echo $1`
ENERGY=`echo $2`
PART=`echo $3`

echo ""
echo ">> `/bin/date` Submitting Condor job(s)..."
#==========================================================
#==========================================================
#==========================================================
echo ">> `/bin/date` PythiaMonoJet: parameters in $1"

mkdir ${DIR} -p
cat ./hepmcanalyzer_cfg_tmp_split.py \
| sed -e s/ENERGYIN/${ENERGY}/ \
| sed -e s/PNUM/${PART}/ \
> ./${DIR}/hepmcanalyzer_cfg_${ENERGY}_split.py

cmsRun ./${DIR}/hepmcanalyzer_cfg_${ENERGY}_split.py