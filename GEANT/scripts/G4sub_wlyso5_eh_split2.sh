#!/bin/bash

for ENERGY in 100 225
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh wlyso5_e wlyso5_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh wlyso5_pi wlyso5_pion_temp_split.in ${ENERGY} ${PART}
    done
  done
