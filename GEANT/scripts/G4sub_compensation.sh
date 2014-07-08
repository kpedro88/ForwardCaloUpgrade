#!/bin/bash

for DET in hcal pbwo4 wlyso
  do
    for PID in elec pion
      do
        for ENERGY in 20 30 50
          do
            ./G4tempE.sh ${DET}_only_${PID} ${DET}_only_${PID}_temp.in ${ENERGY}
          done
      done
  done
  
for DET in hcal pbwo4 wlyso
  do
    for PID in elec pion
      do
        for ENERGY in 100 150 225 300
          do
            for PART in 1 2 3 4 5
              do
                ./G4tempEsplit.sh ${DET}_only_${PID} ${DET}_only_${PID}_temp_split.in ${ENERGY} ${PART}
              done
          done
      done
  done
