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
	temp=np.random.uniform(22.8,24.8,100)
	temp={f:temp for f in filter_set}
	temp['id_galaxy']=np.arange(100)
	
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
		print(f'{ra_region[i]:.2f},{dec_region[i]:.2f}')
		for j in range(len(zrange)):
			temp['ra']=ra_region[i]
			temp['dec']=dec_region[i]
			temp['redshift, cosmological']=(zrange[j][0]+zrange[j][1])/2
			data.append(pd.DataFrame(temp))
	data=pd.concat(data)
	data.to_csv(f'{file_path}test_dataset.csv',index=False)
