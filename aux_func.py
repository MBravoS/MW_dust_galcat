'''Auxiliary functions'''

########################################
# Extinction vector calculation
########################################
def dust_vector(data,band_sel,band_1,band_2,z,A_sel,A_b1,E_b1b2,mag_lim):
	import numpy as np
	import pandas as pd
	
	zsel=(data['zobs']>z[0])&(data['zobs']<z[1])
	base_n=1.0*len(data)
	base_m=np.median(data[band_1])
	base_c=np.median(data[band_1]-data[band_2])
	n_dust,m_dust,c_dust=[[],[],[]]
	for i in range(len(A_sel)):
		mag_sel=(data[band_sel]+A_sel[i])<mag_lim
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
def Pix_Stats_Calc(nside,fname,zcond,band_1,band_2,cc_bands,r_band,r_cut,ebv_map,sigma_cut=False):
	import numpy as np
	import pandas as pd
	import scipy.stats as stats
	
	#No dust
	nside_str=str(int(np.log2(nside))).zfill(2)
	nbin='bin, angular, '+nside_str
	cols=list(set([nbin,band_1,band_2,r_band,'redshift, cosmological']+[b for cc in cc_bands for b in cc]))
	data=pd.read_csv(fname,usecols=cols)
	data=data.loc[zcond]
	bin_id,bin_id_pos,counts=np.unique(data[nbin],return_counts=True,return_inverse=True)
	histogram_bins=np.empty(len(bin_id)+1)
	histogram_bins[:-1]=bin_id-0.5
	histogram_bins[1:]=bin_id+0.5
	mag,temp1,temp2=stats.binned_statistic(data[nbin],data[band_1],statistic='median',bins=histogram_bins)
	colour,temp1,temp2=stats.binned_statistic(data[nbin],data[band_1]-data[band_2],statistic='median',bins=histogram_bins)
	redshift,temp1,temp2=stats.binned_statistic(data[nbin],data['redshift, cosmological'],statistic='median',bins=histogram_bins)
	colour_comp_1,temp1,temp2=stats.binned_statistic(data[nbin],data[cc_bands[0][0]]-data[cc_bands[0][1]],statistic='median',bins=histogram_bins)
	colour_comp_2,temp1,temp2=stats.binned_statistic(data[nbin],data[cc_bands[1][0]]-data[cc_bands[1][1]],statistic='median',bins=histogram_bins)
	del temp1,temp2
	#With dust
	fname_lc_pos=fname.rfind('/')+1
	lc_name=fname[fname_lc_pos:fname_lc_pos+4]
	bands=[lc_name+i+'-band, apparent' for i in [' u',' g',' r',' i',' z']]
	sigma_band=np.array([6.0399,2.0000,1.6635,3.168,6.6226])*1e-12
	sigma={bands[i]:sigma_band[i] for i in xrange(5)}
	EBV_pix=np.random.choice(ebv_map,size=len(bin_id),replace=False)
	mag_A1,col_E=ExtinctionLaw(EBV_pix,band_1,band_2)
	mag_A2=mag_A1-col_E
	mag_R,temp=ExtinctionLaw(EBV_pix,r_band,band_2)
	del EBV_pix,temp
	data[band_1]+=np.array([mag_A1[i] for i in bin_id_pos])
	data[band_2]+=np.array([mag_A2[i] for i in bin_id_pos])
	data[r_band]+=np.array([mag_R[i] for i in bin_id_pos])
	if sigma_cut:
		sigma_cuts=(data[band_1]<=-2.5*np.log10(sigma_cut*sigma[band_1]))
		sigma_cuts&=(data[band_2]<=-2.5*np.log10(sigma_cut*sigma[band_2]))
	else:
		sigma_cuts=(data[band_1]<=98)
		sigma_cuts&=(data[band_2]<=98)
	mag_cut=(sigma_cuts)&(data[r_band]<=r_cut)
	if np.sum(mag_cut)>0:
		data=data.loc[mag_cut]
		bin_id_AE,counts_AE=np.unique(data[nbin],return_counts=True)
		histogram_bins_AE=np.empty(len(bin_id_AE)+1)
		histogram_bins_AE[:-1]=bin_id_AE-0.5
		histogram_bins_AE[1:]=bin_id_AE+0.5
		mag_A,temp1,temp2=stats.binned_statistic(data[nbin],data[band_1],statistic='median',bins=histogram_bins_AE)
		col_E,temp1,temp2=stats.binned_statistic(data[nbin],data[band_1]-data[band_2],statistic='median',bins=histogram_bins_AE)
	else:
		counts_AE,mag_A,col_E=[np.array([])]*3
	del data
	#mag_A1=np.array([mag_A1[i] for i in bin_id_pos])+data[band_1]
	#mag_A2=np.array([mag_A2[i] for i in bin_id_pos])+data[band_2]
	#mag_R=np.array([mag_R[i] for i in bin_id_pos])+data[r_band]
	#if sigma_cut:
	#	sigma_cuts=(mag_A1<=-2.5*np.log10(sigma_cut*sigma[band_1]))
	#	sigma_cuts&=(mag_A2<=-2.5*np.log10(sigma_cut*sigma[band_2]))
	#else:
	#	sigma_cuts=(mag_A1<=98)
	#	sigma_cuts&=(mag_A2<=98)
	#mag_cut=(sigma_cuts)&(mag_R<=r_cut)
	#if np.sum(mag_cut)>0:
	#	temp,counts_AE=np.unique(data[nbin].loc[mag_cut],return_counts=True)
	#	mag_A1=mag_A1[mag_cut]
	#	mag_A2=mag_A2[mag_cut]
	#	mag_A,temp1,temp2=stats.binned_statistic(data[nbin].loc[mag_cut],mag_A1,statistic='median',bins=histogram_bins)
	#	col_E,temp1,temp2=stats.binned_statistic(data[nbin].loc[mag_cut],mag_A1-mag_A2,statistic='median',bins=histogram_bins)
	#else:
	#	counts_AE,mag_A,col_E=[np.array([])]*3
	#del data
	return(counts,mag,colour,redshift,colour_comp_1,colour_comp_2,bin_id,counts_AE,mag_A,col_E)

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
