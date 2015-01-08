#!/bin/bash

./G4tempE.sh hcal_nonuni3_pion hcal_nonuni3_pion_temp.in 50

for ENERGY in 20 50 100
  do
    ./G4tempE.sh hcal_nonuni3_jet hcal_nonuni3_hepmc_tmp.in ${ENERGY}
  done

for ENERGY in 250 500 750
  do
    for PART in 1 2 3 4 5
      do
        ./G4tempEsplit.sh hcal_nonuni3_jet hcal_nonuni3_hepmc_tmp_split.in ${ENERGY} ${PART}
      done
  done
