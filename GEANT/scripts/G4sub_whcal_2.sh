#!/bin/bash

for ENERGY in 150 225 300
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh whcal_only whcal_only_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh wlyso_whcal_e wlyso_whcal_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh wlyso_whcal_pi wlyso_whcal_pion_temp_split.in ${ENERGY} ${PART}
    done
  done

for ENERGY in 20 50 100
  do
    ./G4tempE.sh wlyso_whcal_jet wlyso_whcal_hepmc_temp.in ${ENERGY}
  done

for PART in 1 2 3 4 5
  do
    ./G4tempEsplit.sh wlyso_whcal_jet wlyso_whcal_hepmc_temp_split.in 500 ${PART}
  done
