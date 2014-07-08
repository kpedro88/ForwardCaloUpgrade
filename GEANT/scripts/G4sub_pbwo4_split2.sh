#!/bin/bash

for ENERGY in 100 225
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh hcal_only hcal_only_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh hcal_pbwo4_e hcal_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh hcal_pbwo4_pi hcal_pion_temp_split.in ${ENERGY} ${PART}
    done
  done
