'''Auxiliary functions'''

########################################
# Extinction vector calculation
########################################
def dust_vector(data,band_sel,band_1,band_2,ebv_test,mag_sel_lim,mag_1_lim,mag_2_lim):
	import numpy as np
	import pandas as pd
	
	mag_filt_list=[k for k in data.columns.values if k[-3:]=='_ap_nodust']
	
	mag_sel=data[band_sel]<mag_sel_lim
	mag_sel&=data[band_1]<mag_1_lim
	mag_sel&=data[band_2]<mag_2_lim
	
	data=data.loc[mag_sel]
	
	base_n=1.0*len(data)
	base_m=np.median(data[band_1])
	base_c=np.median(data[band_1]-data[band_2])
	n_dust,m_dust,c_dust=[[],[],[]]
	
	for i in range(len(ebv_test)):
		dusted_data=data.copy()
		for mfl in mag_filt_list:
			A_mfl,temp=extinction_law(ebv_test[i],mfl,band_1)
			dusted_data[mfl]+=A_mfl
			
		mag_sel=(dusted_data[band_sel]<mag_sel_lim)
		mag_sel&=(dusted_data[band_1]<mag_1_lim)
		mag_sel&=(dusted_data[band_2]<mag_2_lim)
		
		new_n=1.0*np.sum(mag_sel)
		new_mag=np.array(dusted_data[band_1])[mag_sel]
		new_col=np.array(dusted_data[band_1]-dusted_data[band_2])[mag_sel]
		
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
			neigh=pf.get_all_neighbours(nside[i],pid)
			for k in neigh:
				if k not in pix_list[i]:
					bcheck[csv_data[res[i]]==pid]=True
		csv_data[f'{res[i]}_border']=bcheck
	csv_data.to_csv(fname,index=False)

########################################
# Add observational errors per file
########################################
def magz_err_perfile(fname):
	import numpy as np
	import pandas as pd
	
	data=pd.read_csv(fname)
	
	####################
	# mag errors
	####################
	bands=['u_ap','g_ap','r_ap','i_ap','z_ap']
	sigma_band=np.array([6.0399,2.0000,1.6635,3.168,6.6226])*1e-12
	sigma={bands[i]:sigma_band[i] for i in range(5)}
	for k in bands:
		temp_mag=-2.5*np.log10(10**(-data[k]/2.5)+np.random.normal(loc=0,scale=sigma[k],size=len(data)))
		data.loc[~np.isnan(temp_mag),k]=temp_mag.loc[~np.isnan(temp_mag)]
		data.loc[np.isnan(temp_mag),k]=99
		temp_mag=-2.5*np.log10(10**(-data[f'{k}_nodust']/2.5)+np.random.normal(loc=0,scale=sigma[k],size=len(data)))
		data.loc[~np.isnan(temp_mag),f'{k}_nodust']=temp_mag.loc[~np.isnan(temp_mag)]
		data.loc[np.isnan(temp_mag),f'{k}_nodust']=99
	
	####################
	# z errors
	####################
	data['zobs_sim']=data['zobs']*1.0
	
	zerr_sigma_low=np.random.normal(loc=0,scale=0.02,size=len(data))
	zerr_sigma_high=np.random.normal(loc=0,scale=np.abs(0.02*(1+0.5*(data['i_ap']-25.3))),size=len(data))
	zerr_i_sel=data['i_ap']>25.3
	zerr_sigma=np.where(zerr_i_sel,zerr_sigma_high,zerr_sigma_low)
	data['zobs']=data['zobs_sim']+(1+data['zobs_sim'])*zerr_sigma
	
	zerr_sigma_low=np.random.normal(loc=0,scale=0.02,size=len(data))
	zerr_sigma_high=np.random.normal(loc=0,scale=np.abs(0.02*(1+0.5*(data['i_ap_nodust']-25.3))),size=len(data))
	zerr_i_sel=data['i_ap_nodust']>25.3
	zerr_sigma=np.where(zerr_i_sel,zerr_sigma_high,zerr_sigma_low)
	data['zobs_nodust']=data['zobs_sim']+(1+data['zobs_sim'])*zerr_sigma
	
	data.to_csv(fname)

########################################
# Assign pixel IDs to galaxies
########################################
def pix_id(fname,nside,sfd_map,seed):
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	csv_data=pd.read_csv(fname)
	col_names=[]
	pix_ids=[]
	for i in range(len(nside)):
		col_name=f'n{np.log2(nside[i]):.0f}'
		col_names.append(col_name)
		pix=pf.ang2pix(nside[i],csv_data['ra'],csv_data['dec'],nest=False,lonlat=True)
		csv_data[col_name]=pix
		uniqpix=np.unique(pix)
		np.random.seed(i*1000+seed)
		sfd_map_sample=np.random.choice(sfd_map[i],len(uniqpix),replace=False)
		sfd_map_sample={uniqpix[i]:sfd_map_sample[i] for i in range(len(uniqpix))}
		pix_ebv=np.zeros(len(pix))
		for pid in uniqpix:
			pix_ebv[pix==pid]=sfd_map_sample[pid]
		csv_data[f'{col_name}_SFDmap']=pix_ebv
		mag_filt_list=[k for k in csv_data.columns.values if k[-3:]=='_ap']
		for mfl in mag_filt_list:
			A_mfl,temp=extinction_law(csv_data[f'{col_name}_SFDmap'],mfl,mag_filt_list[0])
			csv_data[f'{mfl}_nodust']=csv_data[mfl]*1.0
			csv_data[mfl]+=A_mfl
		pix_ids.append(uniqpix)
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
	pixel_name=[]
	
	if bcheck:
		print(f'Will enforce borders for {fname.split("/")[-1]}')
	
	for ns in nside:
		print('Assembling dataframe')
		pixel_df={}
		nside_key=f'n{np.log2(ns):.0f}'
		for z in zr:
			z_key=f'z{np.sum(z)/2:.2f}'.replace('.','')
			####################
			# Data read
			####################
			data_sel=data_full.loc[(data_full['zobs']>z[0])&(data_full['zobs']<z[1])].copy()
			if bcheck:
				n_nonborder=np.sum(~data_sel[f'{nside_key}_border'])
				data_sel=data_sel.loc[~data_sel[f'{nside_key}_border']]
			else:
				n_nonborder=1
			
			if n_nonborder>0:
				####################
				# Apply mag limits
				####################
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
				pixel_df[f'{nside_key}_{z_key}_count']=counts
				pixel_df[f'{nside_key}_{z_key}_mag']=mag
				pixel_df[f'{nside_key}_{z_key}_col']=col
				pname=fname.replace('galaxies',f'pixel_{nside_key}')
				pd.DataFrame(pixel_df).to_csv(pname,index=False)
				pixel_name.append(pname)
			else:
				pixel_name.append(None)
	
	return(pixel_name)

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
		return(np.array(map_data))

########################################
# Fit the dust vector as function of EBV
########################################
def slope(dust_comp):
	import numpy as np
	import scipy.optimize as opti
	
	def err(a,x,y):
		return abs(np.mean(a*x-y))
	
	md=opti.leastsq(err,0,args=(dust_comp['EBV'],dust_comp['delta']))[0][0]
	mm=opti.leastsq(err,0,args=(dust_comp['EBV'],dust_comp['mag']))[0][0]
	mc=opti.leastsq(err,0,args=(dust_comp['EBV'],dust_comp['col']))[0][0]
	
	return(md,mm,mc)

def slope2(dust_comp,ebvmap,dusted=False):
	import numpy as np
	
	if dusted:
		dust_comp['deltaEBV']=dust_comp['EBV']
		temp=dust_comp.loc[dust_comp['deltaEBV']<0.1]
	else:
		dust_comp['deltaEBV']=dust_comp['EBV']-np.median(ebvmap)
		temp=dust_comp.loc[(dust_comp['deltaEBV']>-0.05)&(dust_comp['deltaEBV']<0.05)]
	max_debv=temp['deltaEBV']==np.max(temp['deltaEBV'])
	min_debv=temp['deltaEBV']==np.min(temp['deltaEBV'])
	cen_debv=np.abs(temp['deltaEBV'])==np.min(np.abs(temp['deltaEBV']))
	
	md=(float(temp.loc[max_debv,'delta'])-float(temp.loc[min_debv,'delta']))/0.1
	mm=(float(temp.loc[max_debv,'mag'])-float(temp.loc[min_debv,'mag']))/0.1
	mc=(float(temp.loc[max_debv,'col'])-float(temp.loc[min_debv,'col']))/0.1
	
	dust_comp['delta']-=np.mean(temp.loc[cen_debv,'delta'])
	dust_comp['mag']-=np.mean(temp.loc[cen_debv,'mag'])
	dust_comp['col']-=np.mean(temp.loc[cen_debv,'col'])
	
	#dust_comp['delta']/=md
	#dust_comp['mag']/=mm
	#dust_comp['col']/=mc
	
	return(dust_comp,md,mm,mc)
