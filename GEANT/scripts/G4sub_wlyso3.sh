#!/bin/bash

for ENERGY in 1 2 3 5 9 11 15 20 30 50
  do
    ./G4tempE.sh wlyso3_pion wlyso3_pion_temp.in ${ENERGY}
  done

for ENERGY in 100 150 225 300 1000 3000
  do
    for PART in 1 2 3 4 5
      do 
        ./G4tempEsplit.sh wlyso3_pion wlyso3_pion_temp_split.in ${ENERGY} ${PART}
      done
  done

for ENERGY in 10 30 100 300
  do
    ./G4tempE.sh wlyso3_muon wlyso3_muon_temp.in ${ENERGY}
  done

./G4temp.sh wlyso3_elec wlyso3_elec_50.in
