#!/bin/bash
#PBS -V
#PBS -N MWdust_DESI
#PBS -k eo
#PBS -l nodes=1:ppn=56
#PBS -l walltime=02:00:00
#PBS -q newton

##########################

rm /fast_scratch2/mbravo/MWdust_data/*GALFORM* /fast_scratch2/mbravo/MWdust_plots/*GALFORM*
cd /home/mbravo/MW_dust_galcat

D=desi
P=True
O=/fast_scratch2/mbravo/MWdust_data/
PLT=/fast_scratch2/mbravo/MWdust_plots/
BC=True
SDM=False
M=56
MS=r_ap
N1=64
N2=512

python MW_dust_ext_calc.py -d $D -p $P -o $O -plt $PLT -bc $BC -sdm $SDM -m $M -ms $MS -n $N1 $N2
