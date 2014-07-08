#!/bin/bash

DIR=`echo $1`
ENERGY=`echo $2`

echo ""
echo ">> `/bin/date` PythiaMonoJet: parameters in $1"

mkdir ${DIR} -p
cat ./hepmcanalyzer_cfg_tmp.py \
| sed -e s/ENERGYIN/${ENERGY}/ \
> ./${DIR}/hepmcanalyzer_cfg_${ENERGY}.py

cmsRun ./${DIR}/hepmcanalyzer_cfg_${ENERGY}.py
