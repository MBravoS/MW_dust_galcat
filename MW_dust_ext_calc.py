########################################
# Library imports
########################################
import os
#import aux_func
import argparse
import data_func
#import main_func

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
	parser.add_argument('-f','--files_to_prepare',help='The function from aux_func to copy the data into the file scheme that the code uses.',default=None)
	parser.add_argument('-o','--out_dir',help='The output folder.',default='./')
	opts=parser.parse_args()
	
	#if not opts.input:
	#	parser.error('An input is needed')
	
	####################
	# Prepare data
	####################
	if opts.data_files=='test' and opts.prepare_data:
		getattr(data_func,opts.data_files)(opts.out_dir)
	
	
	

if __name__ == '__main__':
	main()
