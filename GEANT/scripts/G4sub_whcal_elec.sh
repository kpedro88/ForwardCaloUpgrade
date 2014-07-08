#!/bin/bash

for ENERGY in 100 150 225 300
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh whcal_only whcal_only_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh wlyso_whcal_e wlyso_whcal_elec_temp_split.in ${ENERGY} ${PART}
    done
  done
