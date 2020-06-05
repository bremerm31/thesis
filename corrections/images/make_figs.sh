#!/bin/bash

cd ~/lts-mini-app

warp_type=("uniform" "polynomial")
ics=("lake_at_rest" "carrier_greenspan")

WORK="/home/bremerm31"

for wt in ${warp_type[@]}; do
    for ic in ${ics[@]}; do
        cd sisc2020/$ic/$wt/deva
        source deva_sub.sub

        python scripts/proc_counters.py

    done
done