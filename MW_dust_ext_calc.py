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
def main():
	parser=argparse.ArgumentParser()
	parser.add_argument('-d','--data_files',help='The galaxy catalogue to use. By default creates a simple test dataset.',default='test')
	parser.add_argument('-p','--prepare_data',help='The function from aux_func to copy the data into the file scheme that the code uses.',default=False)
	parser.add_argument('-i','--in_dir',help='The input folder.',default='./')
	parser.add_argument('-o','--out_dir',help='The output folder.',default='./')
	parser.add_argument('-plt','--plot_dir',help='The plot folder.',default='./')
	parser.add_argument('-m','--multithread',help='The number of threads for parallelisation.',default=False)
	parser.add_argument('-bc','--border_check',help='Keep or discard border pixels.',default=False)
	parser.add_argument('-sm','--simple_dust_map',help='Test recovery of homogeneous E(B-V)=0.1 map.',default=False)
	opts=parser.parse_args()
	
	#if not opts.input:
	#	parser.error('An input is needed')
	
	####################
	# Prepare data
	####################
	if opts.prepare_data:
		fnames=getattr(data_func,opts.data_files)(opts.out_dir)
	else:
		fnames=glob.glob(f'{opts.in_dir}/galaxies_{opts.data_files}*csv')
	
	####################
	# Dust vector
	####################
	main_func.dust_vector(fnames,'u_ap','u_ap','z_ap',opts.out_dir,opts.plot_dir)
		
	####################
	# Pixelate data
	####################
	main_func.pixel_assign(fnames,nside=[64],border_check=opts.border_check,multithread=opts.multithread,
							simple_ebv=opts.simple_dust_map)
	
	####################
	# Dust vector
	####################
	main_func.dust_vector(fnames,'u_ap','u_ap','z_ap',opts.out_dir,opts.plot_dir,dusted=True)
	
	####################
	# Pixel properties
	####################

if __name__ == '__main__':
	main()
