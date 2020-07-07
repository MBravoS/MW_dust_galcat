#!/bin/bash
#PBS -V
#PBS -N LC_plots
#PBS -k eo
#PBS -l nodes=1:ppn=56
#PBS -l walltime=04:00:00
#PBS -q newton

##########################

cd /home/mbravo/MW_dust_galcat
python lc_description_plots.py
