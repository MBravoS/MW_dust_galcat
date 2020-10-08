'''Here are the core functions for the map analysis'''

########################################
# Dust vector calculation
########################################
def dust_vector(fnames,band_sel,band_1,band_2,data_dir,plot_dir,zrange,mag_cut=24.8,b1_cut=99,b2_cut=99,multithread=False):
	import time
	import aux_func
	import numpy as np
	import pandas as pd
	import splotch as sp
	import matplotlib.cm as cm
	import multiprocessing as mp
	import matplotlib.pyplot as plot
	
	sp.use_style('splotch.style')
	fs=np.array(plot.rcParams.get('figure.figsize'))
	
	####################
	# Reading data in
	####################
	t0=time.time()
	print('Reading galaxy data for dust vector calculation')
	
	run_name=fnames[0].split("galaxies_")[1].split('_')[0]
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(zrange)))
		temp=[pool.apply(aux_func.zsplit,(f,zrange)) for f in fnames]
	else:
		temp=[aux_func.zsplit(f,zrange) for f in fnames]
	data=[[temp[j][i] for j in range(len(fnames))] for i in range(len(zrange))]
	
	####################
	# Extinction effect
	####################
	EBV=np.linspace(0,0.2,num=200)
	A_sel,temp=aux_func.extinction_law(EBV,band_sel,band_1)
	A_b1,E_b1b2=aux_func.extinction_law(EBV,band_1,band_2)
	A_b2=A_b1-E_b1b2
	
	print('Calculating dust vectors')
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(A_sel)))
		vector_comp=[pool.apply(aux_func.dust_vector,(d,f'{band_sel}_nodust',f'{band_1}_nodust',
														f'{band_2}_nodust',EBV,mag_cut,b1_cut,b2_cut,)) for d in data]
	else:
		vector_comp=[aux_func.dust_vector(d,f'{band_sel}_nodust',f'{band_1}_nodust',
											f'{band_2}_nodust',EBV,mag_cut,b1_cut,b2_cut) for d in data]
	
	delta=[]
	mag=[]
	col=[]
	for i in range(len(zrange)):
		delta.append(vector_comp[i]['delta'])
		mag.append(vector_comp[i]['mag'])
		col.append(vector_comp[i]['col'])
		zstr=f'{(zrange[i][0]+zrange[i][1])/2:.2f}'.translate({ord(c): None for c in '.'})
		csv_name=f'{data_dir}dust_vector_{run_name}_z{zstr}.csv'
		vector_comp[i]['EBV']=EBV
		vector_comp[i].to_csv(csv_name,index=False)
	print('Dust vectors saved')
	
	####################
	# Plots
	####################
	#EBV=[EBV]*len(zrange)
	#vir=cm.get_cmap('viridis')
	#cr=vir(np.array([i/5.0 for i in range(len(zrange))]))
	#zlabel=[f'$z_{{{z[0]},{z[1]}}}$' for z in zrange]
	#
	#print('Making E(B-V) vs delta/mag/col plot')
	#
	#plot.figure(figsize=[fs[0],fs[1]*3])
	#
	#plot.subplot(3,1,1)
	#sp.plot(EBV,delta,c=cr,ylabel='$\Delta\log(\delta+1)$',title=run_name,plabel=zlabel)
	#plot.subplot(3,1,2)
	#sp.plot(EBV,mag,c=cr,ylabel=f'$\Delta {band_1[0]}$ [mag]')
	#sp.plot(EBV[0],A_b1,c='k',linestyle='dashed',plabel=f'A$({band_1[0]})$ [mag]')
	#plot.subplot(3,1,3)
	#sp.plot(EBV,col,c=cr,xlabel='E(B-V) [mag]',ylabel=f'$\Delta({band_1[0]}-{band_2[0]})$ [mag]')
	#sp.plot(EBV[0],E_b1b2,c='k',linestyle='dashed',plabel=f'E$({band_1[0]}-{band_2[0]})$ [mag]')
	#plot.tight_layout()
	#plot.savefig(f'{plot_dir}dust_vector_{run_name}.pdf')
	#plot.close()
	
	####################
	# Split delta
	####################
	print('Calculating delta-split dust vectors')
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(A_sel)))
		split_vector_comp=[pool.apply(aux_func.split_dust_vector,(d,f'{band_sel}_nodust',f'{band_1}_nodust',
														f'{band_2}_nodust',EBV,mag_cut,b1_cut,b2_cut,)) for d in data]
	else:
		split_vector_comp=[aux_func.split_dust_vector(d,f'{band_sel}_nodust',f'{band_1}_nodust',
											f'{band_2}_nodust',EBV,mag_cut,b1_cut,b2_cut) for d in data]
	
	for i in range(len(zrange)):
		zstr=f'{(zrange[i][0]+zrange[i][1])/2:.2f}'.translate({ord(c): None for c in '.'})
		csv_name=f'{data_dir}dust_vector_split_{run_name}_z{zstr}.csv'
		split_vector_comp[i]['EBV']=EBV
		split_vector_comp[i].to_csv(csv_name,index=False)
	print('Split-delta dust vectors saved')
	
	t1=time.time()
	print(f'Dust vectors created in {t1-t0} s')
	
	return(vector_comp)

########################################
# Add observational errors to galaxies
########################################
def magz_err(fnames,nside,multithread):
	import aux_func
	import multiprocessing as mp
	
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(fnames)))
		temp=[pool.apply_async(aux_func.magz_err_perfile,(f,nside,)) for f in fnames]
		temp=[t.get() for t in temp]
	else:
		temp=[aux_func.magz_err_perfile(f,nside) for f in fnames]

########################################
# Assign the galaxies to HEALPix pixels
########################################
def pixel_assign(fnames,nside,border_check=False,simple_ebv=True,multithread=False):
	import time
	import aux_func
	import numpy as np
	import healpy as hp
	import pandas as pd
	import multiprocessing as mp
	import healpy.pixelfunc as pf
	
	####################
	# Reading data in
	####################
	t0=time.time()
	print('Reading Schlegel map')
	if simple_ebv:
		ebv_map=[np.linspace(0,0.2,pf.nside2npix(n)) for n in nside]
	else:
		ebv_map=[aux_func.sfd_map(percent=False,resample=n) for n in nside]
	
	print('Reading galaxy data for pixelisation')
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(fnames)))
		temp=[pool.apply_async(aux_func.pix_id,(fnames[j],nside,ebv_map,j,)) for j in range(len(fnames))]
		temp=[t.get() for t in temp]
	else:
		temp=[aux_func.pix_id(fnames[j],nside,ebv_map,j) for j in range(len(fnames))]
	
	res=temp[0][0]
	pix_ids=[t[1] for t in temp]
	pix_ids=[np.concatenate([p[i] for p in pix_ids]) for i in range(len(nside))]
	
	####################
	# Border check
	####################
	if border_check:
		print('Checking borders')
		if multithread:
			pool=mp.Pool(processes=min(multithread,len(fnames)))
			temp=[pool.apply_async(aux_func.find_border,(f,pix_ids,nside,res,)) for f in fnames]
			temp=[t.get() for t in temp]
		else:
			temp=[aux_func.find_border(f,pix_ids,nside,res) for f in fnames]
	
	t1=time.time()
	print(f'Pixelisation finished in {t1-t0} s')

########################################
# Calculate pixel properties
########################################
def pixel_stat(fnames,nside,band_sel,band_1,band_2,zrange,mag_cut=24.8,b1_cut=99,b2_cut=99,border_check=False,
				intrinsic=False,multithread=False):
	import time
	import aux_func
	import numpy as np
	import pandas as pd
	import multiprocessing as mp
	
	####################
	# Stats calculation
	####################
	t0=time.time()
	print('Calculating pixelised properties')
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(fnames)))
		results=[pool.apply_async(aux_func.pix_stat,(f,nside,band_sel,band_1,band_2,mag_cut,b1_cut,b2_cut,
														zrange,border_check,intrinsic,)) for f in fnames]
		results=[r.get() for r in results]
		pool.close()
	else:
		results=[aux_func.pix_stat(f,nside,band_sel,band_1,band_2,mag_cut,b1_cut,b2_cut,zrange,border_check,intrinsic) for f in fnames]
	
	results=[r2 for r1 in results for r2 in r1 if r2 is not None]
	for r in results:
		if ' failed' in r:
			print(r)
	results=[r for r in results if ' failed' not in r]
	
	t1=time.time()
	print(f'Pixel statistics finished in {t1-t0} s')
	return(sorted(results))
