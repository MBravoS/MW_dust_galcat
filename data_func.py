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
def desi(file_path):
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	####################
	# Reading CSV
	####################
	var_names=['id_galaxy_sam','id_galaxy_sky','zcos','zobs','BCDM','ra','dec','u_ap','u_ab','g_ap','g_ab','r_ap','r_ab','i_ap','i_ab','z_ap','z_ab']
	data=pd.read_csv(f'{file_path}DESI.csv',header=None,names=var_names)
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

########################################
# Read the LSST LC made with Buzzard
########################################
def lsst(file_path):
	import numpy as np
	import pandas as pd
	from astropy import table
	from astropy.io import fits
	import healpy.pixelfunc as pf
	
	####################
	# FITS function
	####################
	def fits2csv(fits_name,data_dir,lc)
		hdulist=fits.open(fits_name)
		tbdata=hdulist[1].data
		ra=tbdata['RA']
		dec=tbdata['DEC']
		z=tbdata['Z']
		absmag=tbdata['AMAG']
		appmag=tbdata['OMAG']
		d={'redshift, cosmological' : z.astype('float64'),
		'ra' : ra.astype('float64'),
		'dec' : dec.astype('float64'),
		'LSST u-band, absolute' : absmag[:,0].astype('float64'),
		'LSST g-band, absolute' : absmag[:,1].astype('float64'),
		'LSST r-band, absolute' : absmag[:,2].astype('float64'),
		'LSST i-band, absolute' : absmag[:,3].astype('float64'),
		'LSST z-band, absolute' : absmag[:,4].astype('float64'),
		'LSST u-band, apparent' : appmag[:,0].astype('float64'),
		'LSST g-band, apparent' : appmag[:,1].astype('float64'),
		'LSST r-band, apparent' : appmag[:,2].astype('float64'),
		'LSST i-band, apparent' : appmag[:,3].astype('float64'),
		'LSST z-band, apparent' : appmag[:,4].astype('float64')}
		fit_table=table.Table(d)
		data=fit_table.to_pandas()
		hdulist.close()
		big_pixels=pf.ang2pix(32,data['ra'],data['dec'],nest=False,lonlat=True)
		pixel_list=np.unique(big_pixels)
		
		return(data,zm_cut,sigma_cuts,big_pixels,pixel_list)
	
	
	
	
	####################
	# Reading FITS
	####################
	fits_list=glob.glob(f'{file_path}*.fit')
	pool=mp.Pool(processes=min(56,len(fits_list)))
	csv_out=[pool1.apply(base.fits2csv,(f,data_dir,lc,mag_cut,band_cut,z_cut,bands,[band_1,band_2],sigma,sigma_cut,)) for f in fits_list]
	data=pd.concat([f[0] for f in csv_out],ignore_index=True)
	zm_cut=np.concatenate([f[1] for f in csv_out])
	sigma_cuts=np.concatenate([f[2] for f in csv_out])
	big_pixels=np.concatenate([f[3] for f in csv_out])
	big_pixel_ID=np.unique(np.concatenate([f[4] for f in csv_out]))
