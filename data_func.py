'''Here is where the user should add the functions that they need to
prepare the data on the file style expected by the main function.
The main code expects the catalogue to be in CSV files, splitted by
HEALPix pixelization, each contaning at least the redshift and two
bands for the photometry'''

########################################
# Create a simple test mock
########################################
def test(file_path):
	import os
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	####################
	# DF definition
	####################
	filter_set=['u_ap','u_ab','z_ap','z_ab']
	zrange=[[0.0,0.3],[0.3,0.6],[0.6,0.9],[0.9,1.2],[1.2,2.5]]
	
	####################
	# Create galaxies
	####################
	temp=np.linspace(22.8,24.8,1000)
	temp={f:temp for f in filter_set}
	temp['id_galaxy']=np.arange(1000)
	
	####################
	# Position galaxies
	####################
	'''A ring order is used for this set.'''
	pix_region=pf.get_all_neighbours(32,theta=90,phi=0,nest=False,lonlat=True)
	pix_region=np.concatenate([pix_region,np.array([pf.ang2pix(32,90,0,nest=False,lonlat=True)])])
	pix_map=np.zeros(pf.nside2npix(32))
	pix_map[pix_region]=1
	pix_map=pf.ud_grade(pix_map,64,order_in='RING',order_out='RING')
	pix_region=np.arange(pf.nside2npix(64))[pix_map>0]
	ra_region,dec_region=pf.pix2ang(64,pix_region,nest=False,lonlat=True)
	data=[]
	for i in range(len(pix_region)):
		for j in range(len(zrange)):
			temp['ra']=ra_region[i]
			temp['dec']=dec_region[i]
			temp['zobs']=(zrange[j][0]+zrange[j][1])/2
			data.append(pd.DataFrame(temp))
	data=pd.concat(data)
	big_pixels=pf.ang2pix(32,data['ra'],data['dec'],nest=False,lonlat=True)
	big_pixel_ID=np.unique(big_pixels)
	
	####################
	# File save
	####################
	fnames=[]
	for pix in big_pixel_ID:
		fname=f'{file_path}galaxies_test_{pix}.csv'
		fnames.append(fname)
		data_subset=data.loc[big_pixels==pix]
		data_subset.to_csv(fname,index=False)
	print('Test data created')
	return(fnames)

########################################
# Read the DESI LC made with GALFORM
########################################
def desi(lc,band_1,band_2,test=False,mag_cut=None,band_cut=None,z_cut=None,sigma_cut=False,nested=False):
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	####################
	# Reading CSV
	####################
	var_names=['id_galaxy_sam','id_galaxy_sky','zcos','zobs','BCDM','ra','dec','u_ap','u_ab','g_ap','g_ab','r_ap','r_ab','i_ap','i_ab','z_ap','z_ab']
	data=pd.read_csv('/fast_scratch1/mbravo/DESI/'+lc+'.csv',header=None,names=var_names)
	data=data.drop(columns=['id_galaxy_sam','zcos','BCDM'])
	big_pixels=pf.ang2pix(32,data['ra'],data['dec'],nest=False,lonlat=True)
	big_pixel_ID=np.unique(big_pixels)
	
	####################
	# Saving CSV
	####################
	fnames=[]
	for pix in big_pixel_ID:
		fname=f'{file_path}galaxies_test_{pix}.csv'
		fnames.append(fname)
		data_subset=data.loc[big_pixels==pix]
		data_subset.to_csv(fname,index=False)
	print('DESI/GALFORM data splitted')
	return(fname)
