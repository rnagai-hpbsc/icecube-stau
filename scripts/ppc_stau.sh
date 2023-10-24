#!/bin/bash

for j in `seq 173 $1`; 
do 
    zeroj=`printf '%08d' $j`
    i=`ls data/$2/Stau_prop/Eslope-*_run${zeroj}.i3`
    echo $i; 
    time ./various_stau_ppc_iceinjection.py -o data/$2/Stau_ppc/ --infile $i; 
done; 
