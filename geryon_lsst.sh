#!/bin/bash
#PBS -V
#PBS -N MWdust_LSST
#PBS -k eo
#PBS -l nodes=1:ppn=56
#PBS -l walltime=36:00:00
#PBS -q newton

##########################

rm /fast_scratch2/mbravo/MWdust_data/*Buzzard* /fast_scratch2/mbravo/MWdust_plots/*Buzzard*
cd /home/mbravo/MW_dust_galcat

D=lsst
P=True
O=/fast_scratch2/mbravo/MWdust_data/
PLT=/fast_scratch2/mbravo/MWdust_plots/
BC=True
SDM=False
M=56
MS=r_ap
B1c=27.3
B2c=27.2

python MW_dust_ext_calc.py -d $D -p $P -o $O -plt $PLT -bc $BC -sdm $SDM -m $M -ms $MS -n 64 512 -b1cut $B1c -b2cut $B2c
