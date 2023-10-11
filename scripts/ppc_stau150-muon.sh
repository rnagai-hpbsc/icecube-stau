#!/bin/bash

#for j in `seq 427 $1`; 
#do 
#    zeroj=`printf '%08d' $j`
#    i=`ls data/$2/Stau_prop/Eslope-*_run${zeroj}.i3`
#    echo $i; 
#    time ./timemeas_generate_ppc_staus_iceinjection.py -o data/$2/Stau_ppc/ --infile $i; 
#done; 
for j in `seq 0 $1`; 
do 
    zeroj=`printf '%08d' $j`
    i=`ls data/$2/Muon_prop/Eslope-*_run${zeroj}.i3`
    echo $i; 
    time ./timemeas_generate_ppc_staus_iceinjection.py -o data/$2/Muon_ppc/ --infile $i; 
done
