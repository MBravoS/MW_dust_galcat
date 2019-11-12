'''Here are the core functions for the map analysis'''

########################################
# Data characterization
########################################
def data_char(fnames,band_sel,band_1,band_2,data_dir,z_range=[[0.0,0.3],[0.3,0.6],[0.6,0.9],[0.9,1.2],[1.2,2.5]],
				mag_cut=24.8,sigma_cut=False,multithread=False):
	import aux_func
	import numpy as np
	import pandas as pd
	import splotch as sp
	import astropy.units as u
	import multiprocessing as mp
	from astropy.cosmology import WMAP7 as cosmo
	
	####################
	# Reading data in
	####################
	#colour_index=f'{band_1[0]}{band_2}'
	data=pd.concat([pd.read_csv(f) for f in fnames])
	
	####################
	# Extinction effect
	####################
	EBV=np.linspace(0,0.2,num=200)
	A_sel,temp=aux_func.extinction_law(EBV,band_sel,band_1)
	A_b1,E_b1b2=aux_func.extinction_law(EBV,band_1,band_2)
	A_b2=A_b1-E_b1b2
	
	if multithread:
		cpus=min(multithread,len(Ar))
		pool1=mp.Pool(processes=cpus)
		counts=[np.array([pool1.apply(base.Extinction_Delta,(r_zcut[i],rb1_zcut[i],rb2_zcut[i],Ar[j],Erb1[j],Erb2[j],mag_lim,))
							for j in xrange(len(Ar))]).astype('float64') for i in xrange(len(r_zcut))]
		median_mag=[np.array([pool1.apply(base.Delta_median_mag,(r_zcut[i],rb1_zcut[i],rb2_zcut[i],Ar[j],Erb1[j],Erb2[j],mag_lim,))
								for j in xrange(len(Ar))]).astype('float64') for i in xrange(len(r_zcut))]
		median_col=[np.array([pool1.apply(base.Delta_median_mag,(r_zcut[i],rb1_zcut[i],rb2_zcut[i],Ar[j],Erb1[j],Erb2[j],mag_lim,))
								for j in xrange(len(Ar))]).astype('float64') for i in xrange(len(r_zcut))]
		counts=[np.log10(c/np.max(c)) for c in counts]
	else:
		vector_comp=[aux_func.dust_vector(data,band_sel,band_1,band_2,z,A_sel,A_b1,E_b1b2,mag_cut) for z in z_range]
	for i in range(len(z_range)):
		vector_comp[i].to_csv(f'{data_dir}dust_vector_comp_z_{(z_range[i][0]+z_range[i][1])/2:.1f}.csv')

########################################
# Assign the galaxies to HEALPix pixels
########################################
def pixel_assign(fnames,nside=[64],border_check=False,multithread=False):
	import aux_func
	import numpy as np
	import pandas as pd
	import healpy as hp
	import multiprocessing as mp
	
	####################
	# Reading data in
	####################
	if multithread:
		cpus=min(multithread,len(fnames))
		pool1=mp.Pool(processes=cpus)
		temp=[pool1.apply_async(aux_func.pix_id,(f,nside,)) for f in fnames]
		temp=[t.get() for t in temp]
	else:
		temp=[aux_func.pix_id(f,nside) for f in fnames]
	res=temp[0][0]
	pix_ids=[t[1] for t in temp]
	pix_ids=[np.concatenate([p[i] for p in pix_ids]) for i in range(len(nside))]
	
	####################
	# Border check
	####################
	if border_check:
		if multithread:
			cpus=min(multithread,len(fnames))
			pool2=mp.Pool(processes=cpus)
			temp=[pool2.apply_async(aux_func.find_border,(f,pix_ids,nside,res,)) for f in fnames]
			temp=[t.get() for t in temp]
		else:
			temp=[aux_func.find_border(f,pix_ids,nside,res) for f in fnames]
	
	return(nside)
