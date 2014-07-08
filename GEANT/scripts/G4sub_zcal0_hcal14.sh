#!/bin/bash

for LAY in 0
  do
    ./G4tempEzcal.sh wlyso_zcal_hcal14_pion wlyso_zcal_hcal14_pion_temp.in 50 ${LAY}
    ./G4tempEzcal.sh wlyso_zcal_hcal14_elec wlyso_zcal_hcal14_elec_temp.in 50 ${LAY}
    for ENERGY in 20 50 100
      do
        ./G4tempEzcal.sh wlyso_zcal_hcal14_jet wlyso_zcal_hcal14_hepmc_temp.in ${ENERGY} ${LAY}
      done
    for PART in 1 2 3 4 5
      do
        ./G4tempEzcalsplit.sh wlyso_zcal_hcal14_jet wlyso_zcal_hcal14_hepmc_temp_split.in 500 ${LAY} ${PART}
      done	  
  done
