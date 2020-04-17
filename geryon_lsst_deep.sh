#!/bin/bash
#PBS -V
#PBS -N MWdust_LSST_deep
#PBS -k eo
#PBS -l nodes=1:ppn=56
#PBS -l walltime=72:00:00
#PBS -q newton

##########################

rm /fast_scratch2/mbravo/MWdust_data/*Buzzard* #/fast_scratch2/mbravo/MWdust_plots/*Buzzard*
cd /home/mbravo/MW_dust_galcat

D=lsst
P=True
O=/fast_scratch2/mbravo/MWdust_data_deep/
PLT=/fast_scratch2/mbravo/MWdust_plots_deep/
BC=True
SDM=False
M=56
MS=r_ap
Mc=26.0
B1c=27.3
B2c=27.2

python MW_dust_ext_calc.py -d $D -p $P -o $O -plt $PLT -bc $BC -sdm $SDM -m $M -ms $MS -mcut $Mc -b1cut $B1c -b2cut $B2c -n 64 512 -z 0.0 0.3 0.6 0.9 1.2 2.5 8.0
