D=test
P=True
O=/home/mbravo/Documents/MWdust_testing_data/
#PLT=
BC=False
SM=True
MS=u_ap

rm /home/mbravo/Documents/MWdust_testing_data/*
python MW_dust_ext_calc.py -d $D -p $P -o $O -plt $O -bc $BC -sm $SM -ms $MS
