#!/bin/bash

for ENERGY in 20 50 100
  do
    ./HMCrunE.sh jet ${ENERGY}
  done

./HMCrunEsplit.sh jet 500 5
