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
	import os
	import glob
	import numpy as np
	import pandas as pd
	import mill_base as base
	import multiprocessing as mp
	
	lc_graph=lc+'_'
	if mag_cut is None:
		data_dir='/fast_scratch2/mbravo/MILL_Data_Full/'
	else:
		data_dir='/fast_scratch2/mbravo/MILL_Data_mCut/'
	lc+='_LC'
	if test==True:
		lc+='_test'
	mag_type='apparent'
	if mag_type not in band_1:
		mag_type='absolute'
	bands=[lc[:4]+i+'-band, apparent' for i in [' u',' g',' r',' i',' z']]
	sigma_band=np.array([6.0399,2.0000,1.6635,3.168,6.6226])*1e-12
	sigma={bands[i]:sigma_band[i] for i in xrange(5)}
	if 'DESI' in lc:
	
	
	
	import healpy.pixelfunc as pf
	####################
	# Reading CSV
	####################
	var_names=['id_galaxy_sam','id_galaxy_sky','zcos','zobs','BCDM','ra','dec','u_ap','u_ab','g_ap','g_ab','r_ap','r_ab','i_ap','i_ab','z_ap','z_ab']
	data=pd.read_csv('/fast_scratch1/mbravo/DESI/'+lc+'.csv',header=None,names=var_names)
	data=data.drop(columns=['id_galaxy_sam','zcos','BCDM'])
	#Preliminary data selection
	if mag_cut is None:
		mag_cut=24.8
		band_cut='DESI r-band, apparent'
	if z_cut is None:
		z_cut=2.5
	for i in xrange(5):
		temp_mag=-2.5*np.log10(10**(-data[bands[i]]/2.5)+np.random.normal(loc=0,scale=sigma_band[i],size=len(data)))
		data.loc[~np.isnan(temp_mag),bands[i]]=temp_mag.loc[~np.isnan(temp_mag)]
		data.loc[np.isnan(temp_mag),bands[i]]=99
	sigma_cuts=np.ones(len(data['redshift, cosmological'])).astype('bool')
	if sigma_cut:
		for b in [band_1,band_2]:
			sigma_cuts&=(data[b]<=-2.5*np.log10(sigma_cut*sigma[b]))
	else:
		for b in [band_1,band_2]:
			sigma_cuts&=(data[b]<=98)
	zerr_sigma_low=np.random.normal(loc=0,scale=0.02,size=len(data))
	zerr_sigma_high=np.random.normal(loc=0,scale=np.abs(0.02*(1+0.5*(data['DESI i-band, apparent']-25.3))),size=len(data))
	zerr_i_sel=data['DESI i-band, apparent']>25.3
	zerr_sigma=np.where(zerr_i_sel,zerr_sigma_high,zerr_sigma_low)
	data['redshift, cosmological']=data['redshift, cosmological']+(1+data['redshift, cosmological'])*zerr_sigma
	zm_cut=(data[band_cut]<=mag_cut)&(data['redshift, cosmological']<=z_cut)
	data=data.loc[zm_cut]
	big_pixels=pf.ang2pix(32,data['ra'].loc[sigma_cuts[zm_cut]],data['dec'].loc[sigma_cuts[zm_cut]],nest=nested,lonlat=True)
	big_pixel_ID=np.unique(big_pixels)
	#Parallel dispersion and galaxies lost due to sigma cut
	ebv=base.Map_Analysis()
	bands_short=[base.Band_Read(i) for i in bands]
	colour_disp=pd.DataFrame(np.empty((5,5)),index=bands_short,columns=bands_short)
	gal_lost=pd.DataFrame(np.empty((5,5)),index=bands_short,columns=bands_short)
	for i in xrange(5):
		for j in xrange(5):
			arrow_dx,arrow_dy=base.ExtinctionLaw(ebv,bands_short[i],bands_short[j])
			arrow_l=(arrow_dx[0]**2+arrow_dy[0]**2)**0.5
			arrow_dx=arrow_dx[0]/arrow_l
			arrow_dy=arrow_dy[0]/arrow_l
			if sigma_cut:
				temp_lost=(data[bands[i]]<=-2.5*np.log10(sigma_cut*sigma[bands[i]]))&(data[bands[j]]<=-2.5*np.log10(sigma_cut*sigma[bands[j]]))
			else:
				temp_lost=(data[bands[i]]<=98)&(data[bands[j]]<=98)
			gal_lost.loc[bands_short[i],bands_short[j]]=1.0*np.sum(temp_lost)/np.sum(zm_cut)
			temp_mag=data[bands[i]].loc[temp_lost]
			if i<j:
				temp_col=data[bands[i]].loc[temp_lost]-data[bands[j]].loc[temp_lost]
			else:
				temp_col=data[bands[j]].loc[temp_lost]-data[bands[i]].loc[temp_lost]
			if i!=j:
				colour_disp.loc[bands_short[i],bands_short[j]]=np.std(temp_col*arrow_dy+temp_mag*arrow_dx)
			else:
				colour_disp.loc[bands_short[i],bands_short[j]]=np.std(temp_mag)
	colour_disp.to_csv('/fast_scratch2/mbravo/MILL_Graph/27_'+lc_graph+mag_type[:3]+'.csv')
	gal_lost.to_csv('/fast_scratch2/mbravo/MILL_Graph/28_'+lc_graph+mag_type[:3]+'.csv')
	#Final data selection
	data=data.loc[sigma_cuts[zm_cut]]
	print 'Fraction of galaxies dropped: '+str(100.0*np.sum(~(sigma_cuts&zm_cut))/len(sigma_cuts))+'%'
	#File save
	fname=[]
	for pix in big_pixel_ID:
		temp_fname=data_dir+lc+'_'+str(pix)+'.csv'
		fname.append(temp_fname)
		data_subset=data.loc[big_pixels==pix]
		data_subset.to_csv(temp_fname,index=False)
	return(fname)
