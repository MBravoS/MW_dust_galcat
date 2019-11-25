'''Auxiliary functions'''

########################################
# Extinction vector calculation
########################################
def dust_vector(data,band_sel,band_1,band_2,A_sel,A_b1,A_b2,E_b1b2,mag_sel_lim,mag_1_lim,mag_2_lim,ebv_map):
	import numpy as np
	import pandas as pd
	
	map_sel,map_A1,map_A2=[0,0,0]
	if ebv_map is not None:
		map_A1,map_E12=extinction_law(ebv_map,band_1,band_2)
		map_sel,temp=extinction_law(ebv_map,band_sel,band_2)
		map_A2=map_A1-map_E12
	
	data[band_sel]+=map_sel
	data[band_1]+=map_A1
	data[band_2]+=map_A2
	
	mag_sel=data[band_sel]<mag_sel_lim
	mag_sel&=data[band_1]<mag_1_lim
	mag_sel&=data[band_2]<mag_2_lim
	
	data=data.loc[mag_sel]
	
	base_n=1.0*len(data)
	base_m=np.median(data[band_1])
	base_c=np.median(data[band_1]-data[band_2])
	n_dust,m_dust,c_dust=[[],[],[]]
	
	for i in range(len(A_sel)):
		
		mag_sel=((data[band_sel]+A_sel[i])<mag_sel_lim)
		mag_sel&=((data[band_1]+A_b1[i])<mag_1_lim)
		mag_sel&=((data[band_2]+A_b2[i])<mag_2_lim)
		
		new_n=1.0*np.sum(mag_sel)
		new_mag=np.array(data[band_1]+A_b1[i])[mag_sel]
		new_col=np.array(data[band_1]-data[band_2]+E_b1b2[i])[mag_sel]
		
		n_dust.append(np.log10(new_n/base_n))
		m_dust.append(np.median(new_mag)-base_m)
		c_dust.append(np.median(new_col)-base_c)
	
	n_dust,m_dust,c_dust=[np.array(n_dust),np.array(m_dust),np.array(c_dust)]
	
	return(pd.DataFrame({'delta':n_dust,'mag':m_dust,'col':c_dust}))

########################################
# Convert E(B-V) to A(b1) and E(b1-b2)
########################################
def extinction_law(ebv,band_1,band_2):
	import numpy as np
	
	R={'u_ap':4.37,'g_ap':3.31,'r_ap':2.32,'i_ap':1.72,'z_ap':1.28}
	A1=R[band_1]*ebv
	A2=R[band_2]*ebv
	E=A1-A2
	
	return(A1,E)

########################################
# Check if pixel has neighbours
########################################
def find_border(fname,pix_list,nside,res):
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	csv_data=pd.read_csv(fname)
	for i in range(len(nside)):
		pids=np.unique(csv_data[res[i]])
		bcheck=np.zeros(len(csv_data)).astype('bool')
		for pid in pids:
			neigh=pf.get_all_neighbours(nside,pid)
			for k in neigh:
				if k not in pix_list[i]:
					bcheck[csv_data[res[i]]==pid]=True
		csv_data[f'{res[i]}_border']=bcheck
	csv_data.to_csv(fname,index=False)

########################################
# Assign pixel IDs to galaxies
########################################
def pix_id(fname,nside,sfd_map):
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	csv_data=pd.read_csv(fname)
	col_names=[]
	pix_ids=[]
	for i in range(len(nside)):
		col_name=f'nside_{np.log2(nside[i]):.0f}'
		col_names.append(col_name)
		pix=pf.ang2pix(nside[i],csv_data['ra'],csv_data['dec'],nest=False,lonlat=True)
		csv_data[col_name]=pix
		csv_data[f'{col_name}_SFDmap']=sfd_map[i][pix]
		pix_ids.append(np.unique(pix))
	csv_data.to_csv(fname,index=False)
	return(col_names,pix_ids)
	
########################################
# Calculate the pixel properties
########################################
def pix_stat(fname,nside,bsel,b1,b2,mcut,b1cut,b2cut,zr,bcheck):
	import numpy as np
	import pandas as pd
	import scipy.stats as stats
	
	data_full=pd.read_csv(fname)
	pixel_df={}
	
	for ns in nside:
		nside_key=f'nside_{np.log2(ns):.0f}'
		for z in zr:
			####################
			# Data read
			####################
			data_sel=data_full.loc[(data_full['zobs']>z[0])&(data_full['zobs']<z[1])].copy()
			if bcheck:
				data_sel=data_sel.loc[~data_sel[f'{nside_key}_border']]
			
			####################
			# Add extinction
			####################
			A1,E12=extinction_law(data_sel[f'{nside_key}_SFDmap'],b1,b2)
			Asel,temp=extinction_law(data_sel[f'{nside_key}_SFDmap'],bsel,b2)
			A2=A1-E12
			
			data_sel[bsel]+=Asel
			data_sel[b1]+=A1
			data_sel[b2]+=A2
			
			mag_sel=data_sel[bsel]<mcut
			mag_sel&=data_sel[b1]<b1cut
			mag_sel&=data_sel[b2]<b2cut
			
			data_sel=data_sel.loc[mag_sel]
			
			####################
			# Pixel values
			####################
			bin_id,bin_id_pos,counts=np.unique(data_sel[nside_key],return_counts=True,return_inverse=True)
			histogram_bins=np.array([b-0.5 for b in bin_id]+[bin_id[-1]+0.5])
			
			mag=stats.binned_statistic(data_sel[nside_key],data_sel[b1],statistic='median',bins=histogram_bins)[0]
			col=stats.binned_statistic(data_sel[nside_key],data_sel[b1]-data_sel[b2],statistic='median',bins=histogram_bins)[0]
			ebv=stats.binned_statistic(data_sel[nside_key],data_sel[f'{nside_key}_SFDmap'],statistic='median',bins=histogram_bins)[0]
			
			pixel_df[f'{nside_key}_pixID']=bin_id
			pixel_df[f'{nside_key}_EBV']=ebv
			pixel_df[f'{nside_key}_z{np.sum(z)/2:.2f}_count'.replace('.','')]=counts
			pixel_df[f'{nside_key}_z{np.sum(z)/2:.2f}_mag'.replace('.','')]=mag
			pixel_df[f'{nside_key}_z{np.sum(z)/2:.2f}_col'.replace('.','')]=col
	
	L=0
	for k in pixel_df.keys():
		L=max(L,len(pixel_df[k]))
	for k in pixel_df.keys():
		if len(pixel_df[k])<L:
			missing=L-len(pixel_df[k])
			pixel_df[k]=np.concatenate([pixel_df[k],np.array([None]*missing)])
	pd.DataFrame(pixel_df).to_csv(fname.replace('galaxies','pixel'),index=False)
	return(fname.replace('galaxies','pixel'))

########################################
# Read in the Schlegel map
########################################
def sfd_map(path='./',percent=True,resample=None,lsst_footprint=True):
	import numpy as np
	import healpy as hp
	import astropy.units as u
	from astropy.io import fits
	from astropy.coordinates import Galactic
	from astropy.coordinates import SkyCoord
	
	fname=path+'SFD_map.fits'
	extinction_map=fits.open(fname)
	map_data=extinction_map[1].data['TEMPERATURE']
	extinction_map.close()
	if resample:
		map_data=hp.ud_grade(map_data,nside_out=resample,order_in='NESTED',order_out='RING')
	else:
		map_data=hp.reorder(map_data,n2r=True)
	if lsst_footprint:
		n=hp.npix2nside(len(map_data))
		b,l=hp.pix2ang(n,np.arange(len(map_data)),lonlat=True)
		c=SkyCoord(b*u.degree,l*u.degree,frame=Galactic)
		c=c.icrs
		dec=c.dec.degree
		map_data=map_data[(dec>-70)&(dec<12.5)&(map_data<=0.2)]
	if percent:
		percentiles=np.percentile(map_data,[75,50,25])
		return(percentiles)
	else:
		return(map_data)
