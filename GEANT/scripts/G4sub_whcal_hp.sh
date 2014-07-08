#!/bin/bash

for INPUT in whcal_only_hp wlyso_whcal_e_hp wlyso_whcal_pi_hp wlyso_whcal_jet_hp
  do
    cat ./turn_on_hp.sh \
    > ${INPUT}.sh
    #chmod +x ${INPUT}.sh
  done

for ENERGY in 20 30 50
  do
    ./G4tempE.sh whcal_only_hp whcal_only_elec_temp.in ${ENERGY}
    ./G4tempE.sh wlyso_whcal_e_hp wlyso_whcal_elec_temp.in ${ENERGY}
    ./G4tempE.sh wlyso_whcal_pi_hp wlyso_whcal_pion_temp.in ${ENERGY}
  done

for ENERGY in 100 150 225 300
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh whcal_only_hp whcal_only_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh wlyso_whcal_e_hp wlyso_whcal_elec_temp_split.in ${ENERGY} ${PART}
        ./G4tempEsplit.sh wlyso_whcal_pi_hp wlyso_whcal_pion_temp_split.in ${ENERGY} ${PART}
    done
  done

for ENERGY in 20 50 100
  do
    ./G4tempE.sh wlyso_whcal_jet_hp wlyso_whcal_hepmc_temp.in ${ENERGY}
  done

for PART in 1 2 3 4 5
  do
    ./G4tempEsplit.sh wlyso_whcal_jet_hp wlyso_whcal_hepmc_temp_split.in 500 ${PART}
  done
