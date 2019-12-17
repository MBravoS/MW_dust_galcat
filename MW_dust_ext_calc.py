########################################
# Library imports
########################################
import os
import glob
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
	parser=argparse.ArgumentParser()
	parser.add_argument('-d','--data_files',help='The galaxy catalogue to use. By default creates a simple test dataset.',default='test')
	parser.add_argument('-p','--prepare_data',help='The function from aux_func to copy the data into the file scheme that the code uses.',
						default=False,type=s2b)
	parser.add_argument('-i','--in_dir',help='The input folder.',default='./')
	parser.add_argument('-o','--out_dir',help='The output folder.',default='./')
	parser.add_argument('-plt','--plot_dir',help='The plot folder.',default='./')
	parser.add_argument('-m','--multithread',help='The number of threads for parallelisation.',default=None,type=int)
	parser.add_argument('-bc','--border_check',help='Keep or discard border pixels.',default=False,type=s2b)
	parser.add_argument('-sm','--simple_dust_map',help='Test recovery of homogeneous E(B-V)=0.1 map.',default=False,type=s2b)
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
	
	#if not opts.input:
	#	parser.error('An input is needed')
	
	####################
	# Prepare data
	####################
	path_dict={'test':'test','desi':'/fast_scratch1/mbravo/DESI/','lsst':'/fast_scratch1/mbravo/LSST/'}
	if opts.prepare_data:
		fnames=getattr(data_func,opts.data_files)(path_dict[opts.data_files],opts.out_dir)
	else:
		fnames=glob.glob(f'{opts.in_dir}/galaxies_{opts.data_files}*csv')
	
	####################
	# Dust vector
	####################
	dust_vector=main_func.dust_vector(fnames,opts.sel_band,opts.band1,opts.band2,opts.out_dir,opts.plot_dir,opts.zbins,
										mag_cut=opts.sel_mag_cut,b1_cut=opts.band1_mag_cut,b2_cut=opts.band2_mag_cut,
										dusted=False,multithread=opts.multithread)
	
	####################
	# Pixelate data
	####################
	main_func.pixel_assign(fnames,opts.nside,border_check=opts.border_check,multithread=opts.multithread,
							simple_ebv=opts.simple_dust_map)
	
	####################
	# Dust vector
	####################
	temp=main_func.dust_vector(fnames,opts.sel_band,opts.band1,opts.band2,opts.out_dir,opts.plot_dir,opts.zbins,
								mag_cut=opts.sel_mag_cut,b1_cut=opts.band1_mag_cut,b2_cut=opts.band2_mag_cut,
								dusted=True,multithread=opts.multithread)
	
	####################
	# Pixel properties
	####################
	pnames=main_func.pixel_stat(fnames,opts.nside,opts.sel_band,opts.band1,opts.band2,opts.zbins,
								mag_cut=opts.sel_mag_cut,b1_cut=opts.band1_mag_cut,b2_cut=opts.band2_mag_cut,
								border_check=opts.border_check,multithread=opts.multithread)
	
	####################
	# Dust map
	####################
	#main_func.dust_mapping(pnames,dust_vector,opts.nside,opts.zbins,opts.out_dir,opts.plot_dir)
	#main_func.dust_mapping2(pnames,dust_vector,opts.nside,opts.zbins,opts.out_dir,opts.plot_dir)
	main_func.dust_mapping3(pnames,dust_vector,opts.nside,opts.zbins,opts.out_dir,opts.plot_dir)

if __name__ == '__main__':
	main()
