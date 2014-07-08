#!/bin/bash

./G4temp.sh wlyso6_elec wlyso6_elec_50.in
./G4tempE.sh wlyso6_pion wlyso6_pion_temp.in 50

for ENERGY in 20 50 100
  do
    ./G4tempE.sh wlyso6_jet hepmc06_tmp.in ${ENERGY}
  done

for PART in 1 2 3 4 5
  do
    ./G4tempEsplit.sh wlyso6_jet hepmc06_tmp_split.in 500 ${PART}
  done
