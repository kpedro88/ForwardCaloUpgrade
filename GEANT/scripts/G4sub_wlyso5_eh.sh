#!/bin/bash

#for ENERGY in 20 30 50
for ENERGY in 20 30
  do
    ./G4tempE.sh wlyso5_e wlyso5_elec_temp.in ${ENERGY}
    ./G4tempE.sh wlyso5_pi wlyso5_pion_temp.in ${ENERGY}
  done

for ENERGY in 150 300
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh wlyso5_e wlyso5_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh wlyso5_pi wlyso5_pion_temp_split.in ${ENERGY} ${PART}
    done
  done
