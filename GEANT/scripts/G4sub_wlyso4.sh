#!/bin/bash

for ENERGY in 11 50
  do
    ./G4tempE.sh wlyso4_pion wlyso4_pion_temp.in ${ENERGY}
  done

./G4temp.sh wlyso4_elec wlyso4_elec_50.in
