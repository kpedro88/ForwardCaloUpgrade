#!/bin/bash

./G4temp.sh wlyso5_elec wlyso5_elec_50.in
./G4tempE.sh wlyso5_pion wlyso5_pion_temp.in 50

for ENERGY in 20 50 100
  do
    ./G4tempE.sh wlyso5_jet hepmc05_tmp.in ${ENERGY}
  done

for PART in 1 2 3 4 5
  do
    ./G4tempEsplit.sh wlyso5_jet hepmc05_tmp_split.in 500 ${PART}
  done
