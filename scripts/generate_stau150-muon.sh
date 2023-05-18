#!/bin/bash

for i in `seq 0 $1`; 
do 
    echo $i; 
    time ./tmp_only_generation-genereate_ppc_staus_iceinjection.py -o data/$2/Stau_prop/Eslope-Stau.i3 -w --minloge 5 --maxloge 12 -r $i -n 200; 
    zeroi=`printf '%08d' $i`
    time ./muon_generate_from_stau_iceinjection.py -o data/$2/Muon_prop/Eslope-Muon.i3 --infile data/$2/Stau_prop/Eslope*_run${zeroi}.i3 -w --minloge 5 --maxloge 12 -r $i; 
done
