########################################
# Library imports
########################################
import os
import glob
import time
import argparse
import data_func
import main_func

########################################
# User variables
########################################

########################################
# Main code
########################################
def s2b(s):
	import distutils.util as u
	
	b=bool(u.strtobool(s))
	return(b)

def main():
	t0=time.time()
	
	parser=argparse.ArgumentParser()
	parser.add_argument('-d','--data_files',help='The galaxy catalogue to use. By default creates a simple test dataset.',default='test')
	parser.add_argument('-i','--in_dir',help='The input folder.',default='./')
	parser.add_argument('-o','--out_dir',help='The output folder.',default='./')
	parser.add_argument('-plt','--plot_dir',help='The plot folder.',default='./')
	parser.add_argument('-m','--multithread',help='The number of threads for parallelisation.',default=None,type=int)
	parser.add_argument('-bc','--border_check',help='Keep or discard border pixels.',default=False,type=s2b)
	parser.add_argument('-sdm','--simple_dust_map',help='Test recovery of homogeneous E(B-V)=0.1 map.',default=False,type=s2b)
	parser.add_argument('-n','--nside',help='HEALPix nside to be used for the map recovery.',default=[64],nargs='+',type=int)
	parser.add_argument('-z','--zbins',help='Bin edges for the redshift bins.',default=None,nargs='+',type=float)
	parser.add_argument('-ms','--sel_band',help='Selection band.',default='r_ap')
	parser.add_argument('-b1','--band1',help='Band 1.',default='u_ap')
	parser.add_argument('-b2','--band2',help='Band 2.',default='z_ap')
	parser.add_argument('-mcut','--sel_mag_cut',help='Magnitude cut to apply to selection band.',default=24.8,type=float)
	parser.add_argument('-b1cut','--band1_mag_cut',help='Magnitude cut to apply to band 1.',default=99,type=float)
	parser.add_argument('-b2cut','--band2_mag_cut',help='Magnitude cut to apply to band 2.',default=99,type=float)
	opts=parser.parse_args()
	
	####################
	# Check inputs
	####################
	if opts.zbins is None:
		opts.zbins=[[0.0,0.3],[0.3,0.6],[0.6,0.9],[0.9,1.2],[1.2,2.5]]
	else:
		opts.zbins=[[opts.zbins[i],opts.zbins[i+1]] for i in range(len(opts.zbins)-1)]
	
	####################
	# Read data
	####################
	path_dict={'test':'test','desi':'/fast_scratch1/mbravo/DESI/','lsst':'/fast_scratch1/mbravo/LSST/'}
	fname_dict={'desi':'GALFORM','lsst':'Buzzard'}
	fnames=glob.glob(f'{opts.in_dir}/galaxies_{fname_dict[opts.data_files]}*csv')
	data_exist=True
	if len(fnames)==0:
		fnames=getattr(data_func,opts.data_files)(path_dict[opts.data_files],opts.out_dir)
		data_exist=False
	
	####################
	# Pixelate data
	####################
	if not data_exist:
		main_func.pixel_assign(fnames,opts.nside,border_check=opts.border_check,multithread=opts.multithread,
								simple_ebv=opts.simple_dust_map)
	
	####################
	# Intrinsic pixels
	####################
	pnames=main_func.pixel_stat(fnames,opts.nside,opts.sel_band,opts.band1,opts.band2,opts.zbins,
								mag_cut=opts.sel_mag_cut,b1_cut=opts.band1_mag_cut,b2_cut=opts.band2_mag_cut,
								border_check=opts.border_check,intrinsic=True,multithread=opts.multithread)
	
	####################
	# Add errors
	####################
	if not data_exist:
		main_func.magz_err(fnames,multithread=opts.multithread)
	
	####################
	# Dust vector
	####################
	dust_vector=glob.glob(f'{opts.in_dir}/dust_vector_{fname_dict[opts.data_files]}*csv')
	if len(dust_vector)>0:
		dust_vector=[pd.read_csv(dv) for dv in sorted(dust_vector)]
	else:
		dust_vector=main_func.dust_vector(fnames,opts.sel_band,opts.band1,opts.band2,opts.out_dir,opts.plot_dir,opts.zbins,
											mag_cut=opts.sel_mag_cut,b1_cut=opts.band1_mag_cut,b2_cut=opts.band2_mag_cut,
											multithread=opts.multithread)
	
	####################
	# Pixel properties
	####################
	pnames=main_func.pixel_stat(fnames,opts.nside,opts.sel_band,opts.band1,opts.band2,opts.zbins,
								mag_cut=opts.sel_mag_cut,b1_cut=opts.band1_mag_cut,b2_cut=opts.band2_mag_cut,
								border_check=opts.border_check,multithread=opts.multithread)
	
	####################
	# Dust map
	####################
	main_func.dust_mapping(pnames,dust_vector,opts.nside,opts.zbins,opts.out_dir,opts.plot_dir)
	
	t1=time.time()
	print(f'Total running time = {t1-t0} s')

if __name__ == '__main__':
	main()
