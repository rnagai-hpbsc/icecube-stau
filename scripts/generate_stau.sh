#!/bin/bash

time for i in `seq 0 $1`; 
do 
    echo $i; 
    time ./various_stau_generate_iceinjection.py -o data/$2/Stau_prop/Eslope-Stau.i3 -w --minloge 5 --maxloge 12 -r $i -n 200 --mass $3; 
    zeroi=`printf '%08d' $i`
    time ./muon_generate_from_stau_iceinjection.py -o data/$2/Muon_prop/Eslope-Muon.i3 --infile data/$2/Stau_prop/Eslope*_run${zeroi}.i3 -w --minloge 5 --maxloge 12 -r $i --mass $3; 
done
