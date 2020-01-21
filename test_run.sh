D=test
P=True
O=/home/mbravo/Documents/MWdust_testing_data/
#PLT=
BC=False
SDM=True
MS=r_ap

rm /home/mbravo/Documents/MWdust_testing_data/*
python MW_dust_ext_calc.py -d $D -p $P -o $O -plt $O -bc $BC -sdm $SDM -ms $MS
