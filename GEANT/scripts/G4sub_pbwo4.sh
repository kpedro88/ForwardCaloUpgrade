#!/bin/bash

for ENERGY in 20 30 50
  do
    ./G4tempE.sh hcal_only hcal_only_elec_temp.in ${ENERGY}
    ./G4tempE.sh hcal_pbwo4_e hcal_elec_temp.in ${ENERGY}
    ./G4tempE.sh hcal_pbwo4_pi hcal_pion_temp.in ${ENERGY}
  done

for ENERGY in 150 300
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh hcal_only hcal_only_elec_temp.in ${ENERGY} ${PNUMBER}
        ./G4tempEsplit.sh hcal_pbwo4_e hcal_elec_temp.in ${ENERGY} ${PNUMBER}
        ./G4tempEsplit.sh hcal_pbwo4_pi hcal_pion_temp.in ${ENERGY} ${PNUMBER}
    done
  done
