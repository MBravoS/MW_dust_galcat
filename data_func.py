'''Here is where the user should add the functions that they need to
prepare the data on the file style expected by the main function.
The main code expects the catalogue to be in CSV files, splitted by
HEALPix pixelization, each contaning at least the redshift and two
bands for the photometry'''

########################################
# Create a simple test mock
########################################
def test(temp,file_path):
	import time
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	####################
	# DF definition
	####################
	filter_set=['u','g','r','i','z']
	filter_set=[f'{fs}_ap' for fs in filter_set]+[f'{fs}_ab' for fs in filter_set]
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
def desi(in_path,out_path):
	import time
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	####################
	# Reading CSV
	####################
	t0=time.time()
	var_names=['id_galaxy_sam','id_galaxy_sky','zcos','zobs','BCDM','ra','dec','u_ap','u_ab','g_ap','g_ab','r_ap','r_ab','i_ap','i_ab','z_ap','z_ab']
	data=pd.read_csv(f'{in_path}DESI_LC.csv',header=None,names=var_names)
	data=data.drop(columns=['id_galaxy_sam','zcos','BCDM'])
	big_pixels=pf.ang2pix(32,data['ra'],data['dec'],nest=False,lonlat=True)
	big_pixel_ID=np.unique(big_pixels)
	
	####################
	# Saving CSV
	####################
	fnames=[]
	for pix in big_pixel_ID:
		fname=f'{out_path}galaxies_GALFORM_{pix}.csv'
		fnames.append(fname)
		data_subset=data.loc[big_pixels==pix]
		data_subset.to_csv(fname,index=False)
	
	t1=time.time()
	print(f'Data prepared in {t1-t0} s')
	return(fnames)

########################################
# Read the LSST LC made with Buzzard
########################################

####################
# FITS function
####################
def fits2csv(fits_name):
	import numpy as np
	import pandas as pd
	from astropy import table
	from astropy.io import fits
	import healpy.pixelfunc as pf
	
	hdulist=fits.open(fits_name)
	tbdata=hdulist[1].data
	ra=tbdata['RA']
	dec=tbdata['DEC']
	z=tbdata['Z']
	absmag=tbdata['AMAG']
	appmag=tbdata['OMAG']
	d={'zobs':z.astype('float64'),
		'ra':ra.astype('float64'),
		'dec':dec.astype('float64'),
		'u_ab':absmag[:,0].astype('float64'),
		'g_ab':absmag[:,1].astype('float64'),
		'r_ab':absmag[:,2].astype('float64'),
		'i_ab':absmag[:,3].astype('float64'),
		'z_ab':absmag[:,4].astype('float64'),
		'u_ap':appmag[:,0].astype('float64'),
		'g_ap':appmag[:,1].astype('float64'),
		'r_ap':appmag[:,2].astype('float64'),
		'i_ap':appmag[:,3].astype('float64'),
		'z_ap':appmag[:,4].astype('float64')}
	fit_table=table.Table(d)
	data=fit_table.to_pandas()
	hdulist.close()
	data['big_pix']=pf.ang2pix(32,theta=data['ra'],phi=data['dec'],lonlat=True,nest=False)
	return(data)

####################
# Read function
####################

def lsst(in_path,out_path):
	import glob
	import time
	import numpy as np
	import pandas as pd
	import multiprocessing as mp
	
	####################
	# Reading FITS
	####################
	t0=time.time()
	
	fits_list=glob.glob(f'{in_path}*.fit')
	pool=mp.Pool(processes=min(56,len(fits_list)))
	df_out=[pool.apply_async(fits2csv,(f,)) for f in fits_list]
	df_out=pd.concat([c.get() for c in df_out])
	
	df_names=[]
	bigpix=np.unique(df_out['big_pix'])
	for bp in bigpix:
		df_temp=df_out.loc[df_out['big_pix']==bp]
		if len(df_temp)>100:
			df_temp.drop(columns=['big_pix'],inplace=True)
			df_name=f'{out_path}galaxies_Buzzard_{bp}.csv'
			df_temp.to_csv(df_name,index=False)
			df_names.append(df_name)
	
	t1=time.time()
	print(f'Data prepared in {t1-t0} s')
	
	return(df_names)
