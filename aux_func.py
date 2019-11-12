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
def pix_id(fname,nside):
	import numpy as np
	import pandas as pd
	import healpy.pixelfunc as pf
	
	csv_data=pd.read_csv(fname)
	col_names=[]
	pix_ids=[]
	for n in nside:
		col_name=f'nside_{np.log2(n):.0f}'
		col_names.append(col_name)
		pix=pf.ang2pix(n,csv_data['ra'],csv_data['dec'],nest=False,lonlat=True)
		csv_data[col_name]=pix
		pix_ids.append(np.unique(pix))
	csv_data.to_csv(fname,index=False)
	return(col_names,pix_ids)
