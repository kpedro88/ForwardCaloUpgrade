#!/bin/bash

./G4temp.sh wlyso7_elec wlyso7_elec_50.in
./G4tempE.sh wlyso7_pion wlyso7_pion_temp.in 50

for ENERGY in 20 50 100
  do
    ./G4tempE.sh wlyso7_jet hepmc07_tmp.in ${ENERGY}
  done

for PART in 1 2 3 4 5
  do
    ./G4tempEsplit.sh wlyso7_jet hepmc07_tmp_split.in 500 ${PART}
  done
