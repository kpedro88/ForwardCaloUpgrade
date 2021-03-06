#!/bin/bash

for ENERGY in 1 2 3 5 9 11 15 20 30 50
  do
    ./G4tempE.sh wlyso2_pion wlyso2_pion_temp.in ${ENERGY}
  done

for ENERGY in 100 150 225 300 1000 3000
  do
    for PART in 1 2 3 4 5
      do 
        ./G4tempEsplit.sh wlyso2_pion wlyso2_pion_temp_split.in ${ENERGY} ${PART}
      done
  done

for ENERGY in 10 30 100 300
  do
    ./G4tempE.sh wlyso2_muon wlyso2_muon_temp.in ${ENERGY}
  done

./G4temp.sh wlyso2_elec wlyso2_elec_50.in
