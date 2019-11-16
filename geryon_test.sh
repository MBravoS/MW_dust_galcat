#!/bin/bash
#PBS -V
#PBS -N MWdust_test
#PBS -k eo
#PBS -l nodes=1:ppn=56
#PBS -l walltime=02:00:00
#PBS -q newton
#PBS -m ae
#PBS -M matias.bravo@icrar.org

##########################

#python MW_dust_ext_calc.py -d test -p True -o /home/mbravo/Documents/MWdust_testing_data/ -plt /home/mbravo/Documents/MWdust_testing_data/ -bc True
python MW_dust_ext_calc.py -d test -i /home/mbravo/Documents/MWdust_testing_data/ -o /home/mbravo/Documents/MWdust_testing_data/ -plt /home/mbravo/Documents/MWdust_testing_data/ -bc True
