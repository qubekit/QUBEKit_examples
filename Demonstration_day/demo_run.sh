#!/bin/bash

# Run QUBEKit in an automated way
# first run the qubekit command

QUBEKit -i $1 -config demo.ini -end lennard_jones

# When done move into the finallise folder

IFS='.' read pre suf <<< $1

filename=QUBEKit_${pre}_2019_07_05_JH

cd $filename/11_finalise

# now run the python script to get the density

python ../../../scripts/run_all.py liq $2 -


