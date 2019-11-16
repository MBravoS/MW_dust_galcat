#!/bin/bash
#PBS -V
#PBS -N MWdust_test
#PBS -k eo
#PBS -l nodes=1:ppn=56
#PBS -l walltime=02:00:00
#PBS -q newton
#PBS -e /home/mbravo/MW_dust_galcat/test_error.txt
#PBS -o /home/mbravo/MW_dust_galcat/test_output.txt

##########################

cd /home/mbravo/MW_dust_galcat
which python
python MW_dust_ext_calc.py -d test -p True -o /fast_scratch2/mbravo/MWdust_data/ -plt /fast_scratch2/mbravo/MWdust_plots/ -bc True
