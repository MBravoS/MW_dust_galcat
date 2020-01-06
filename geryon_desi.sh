#!/bin/bash
#PBS -V
#PBS -N MWdust_DESI
#PBS -k eo
#PBS -l nodes=1:ppn=56
#PBS -l walltime=02:00:00
#PBS -q newton

##########################

rm /fast_scratch2/mbravo/MWdust_data/*desi* /fast_scratch2/mbravo/MWdust_plots/*desi*
cd /home/mbravo/MW_dust_galcat

D=desi
P=True
O=/fast_scratch2/mbravo/MWdust_data/
PLT=/fast_scratch2/mbravo/MWdust_plots/
BC=True
SM=True
M=56
MS=r_ap

python MW_dust_ext_calc.py -d $D -p $P -o $O -plt $PLT -bc $BC -sm $SM -m $M -ms $MS
