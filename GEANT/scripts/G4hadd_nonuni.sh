#!/bin/bash

DIR=/data/users/pedrok/FullSim

for ENERGY in 250 500 750
  do
    for FILE in pbwo4/hcal_nonuni_monojet
      do
        hadd ${DIR}/${FILE}_${ENERGY}gev_10k.root ${DIR}/${FILE}_${ENERGY}gev_part*_10k.root
        rm ${DIR}/${FILE}_${ENERGY}gev_part*_10k.root
      done
  done
