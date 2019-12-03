'''Here are the core functions for the map analysis'''
########################################
# Dust map recovery
########################################
def dust_mapping(pnames,dvec,nside,zrange,out_dir,plot_dir):
	import aux_func
	import numpy as np
	import pandas as pd
	import splotch as sp
	import matplotlib.pyplot as plot
	
	for i in range(len(nside)):
		ns=nside[i]
		nside_key=f'n{np.log2(ns):.0f}'
		data=pd.concat([pd.read_csv(p) for p in pnames if nside_key in p],ignore_index=True)
		EBV_recovery=[]
		z_label=[]
		for j in range(len(zrange)):
			zr=zrange[j]
			z_key=f'z{np.sum(zr)/2:.2f}'.replace('.','')
			
			####################
			# Delta calc
			####################
			delta=data[f'{nside_key}_{z_key}_count']+np.random.uniform(-0.5,0.5,len(data))
			data[f'{nside_key}_{z_key}_delta']=np.log10(delta/np.mean(delta))
			
			####################
			# Read dust vector
			####################
			dust_vector=pd.read_csv(f'{out_dir}dust_vector_{pnames[0].split("/")[-1].split("_")[2]}_{z_key}_dusted.csv')
			ebv=dust_vector['EBV']
			ebv2d,ebv2m,ebv2c=aux_func.slope(dust_vector)
			ebv_from_delta=(data[f'{nside_key}_{z_key}_delta']-np.max(data[f'{nside_key}_{z_key}_delta']))/ebv2d
			ebv_from_mag=(data[f'{nside_key}_{z_key}_mag']-np.min(data[f'{nside_key}_{z_key}_mag']))/ebv2m
			ebv_from_col=(data[f'{nside_key}_{z_key}_col']-np.min(data[f'{nside_key}_{z_key}_col']))/ebv2c
			ebv_recover=np.zeros(len(ebv_from_delta))
			for k in range(len(ebv_from_delta)):
				ebv_recover[k]=np.dot([ebv_from_delta[k],ebv_from_mag[k],ebv_from_col[k]],[1,1,1])/(3**0.5)
			#ebv_from_delta=(data[f'{nside_key}_{z_key}_delta']-np.max(data[f'{nside_key}_{z_key}_delta']))
			#ebv_from_mag=(data[f'{nside_key}_{z_key}_mag']-np.min(data[f'{nside_key}_{z_key}_mag']))
			#ebv_from_col=(data[f'{nside_key}_{z_key}_col']-np.min(data[f'{nside_key}_{z_key}_col']))
			#ebv_recover=np.zeros(len(ebv_from_delta))
			#for k in range(len(ebv_from_delta)):
			#	ebv_recover[k]=np.dot([ebv_from_delta[k],ebv_from_mag[k],ebv_from_col[k]],[ebv2d,ebv2m,ebv2c])/(ebv2d**2+ebv2m**2+ebv2c**2)**0.5
			#print((ebv2d**2+ebv2m**2+ebv2c**2)**0.5)
			#ebv_recover/=3.3
			
			#ebv_recover+=np.min(data[f'{nside_key}_EBV'])-np.min(ebv_recover)
			EBV_recovery.append(ebv_recover)
			z_label.append(z_key)
		
		####################
		# Plots
		####################
		plot.figure()
		
		for j in range(len(zrange)):
			sp.scatter(data[f'{nside_key}_EBV'],EBV_recovery[j]+np.min(data[f'{nside_key}_EBV']),plabel=z_label[j],c=f'C{j}',
						xlabel='$E(B-V)_\mathrm{input}$',ylabel='$E(B-V)_\mathrm{recovered}$')
		sp.axline(m=1,plabel='1:1',color='k',linestyle='dashed')
		plot.savefig(f'{plot_dir}ebv_recovery_zbin.pdf')
		plot.figure()
		
		for j in range(len(zrange)):
			sp.scatter(data[f'{nside_key}_EBV']-np.min(data[f'{nside_key}_EBV']),EBV_recovery[j],plabel=z_label[j],c=f'C{j}',
						xlabel='$\Delta E(B-V)_\mathrm{input}$',ylabel='$\Delta E(B-V)_\mathrm{recovered}$')
		sp.axline(m=1,plabel='1:1',color='k',linestyle='dashed')
		plot.savefig(f'{plot_dir}delta_ebv_recovery_zbin.pdf')

########################################
# Dust vector calculation
########################################
def dust_vector(fnames,band_sel,band_1,band_2,data_dir,plot_dir,zrange,mag_cut=24.8,b1_cut=99,b2_cut=99,dusted=False,multithread=False):
	import aux_func
	import numpy as np
	import pandas as pd
	import splotch as sp
	import matplotlib.cm as cm
	import multiprocessing as mp
	import matplotlib.pyplot as plot
	
	sp.use_style('splotch.style')
	fs=np.array(plot.rcParams.get('figure.figsize'))
	fd=''
	
	####################
	# Reading data in
	####################
	print('Reading galaxy data for dust vector calculation')
	run_name=fnames[0].split("galaxies_")[1].split('_')[0]
	data=pd.concat([pd.read_csv(f) for f in fnames])
	ebv=np.zeros(len(data))
	if dusted:
		fd+='_dusted'
		sfd_key=[k for k in data.columns.values if 'SFD' in k]
		if len(sfd_key)>1:
			sfd_nside=np.array([int(s.split('_')[0][1:]) for s in sfd_key])
			maxfound=False
			j=0
			while not maxfound:
				if sfd_nside[j]==np.max(sfd_nside):
					maxfound=True
				else:
					j+=1
			sfd_key=sfd_key[np.arange(len(sfd_key))[sfd_nside==np.max(sfd_nside)][0]]
			ebv=data[sfd_key]
		else:
			ebv=np.array(data[sfd_key[0]])
	
	####################
	# Extinction effect
	####################
	EBV=np.linspace(0,0.2,num=200)
	A_sel,temp=aux_func.extinction_law(EBV,band_sel,band_1)
	A_b1,E_b1b2=aux_func.extinction_law(EBV,band_1,band_2)
	A_b2=A_b1-E_b1b2
	
	if dusted:
		print('Calculating observed dust vectors')
	else:
		print('Calculating intrinsic dust vectors')
	if multithread:
		pool=mp.Pool(processes==min(multithread,len(Ar)))
		vector_comp=[pool.apply(aux_func.dust_vector(data.loc[(data['zobs']>z[0])&(data['zobs']<z[1])].copy(),band_sel,band_1,band_2,
														A_sel,A_b1,A_b2,E_b1b2,mag_cut,b1_cut,b2_cut,
														ebv[(data['zobs']>z[0])&(data['zobs']<z[1])],)) for z in zrange]
	else:
		vector_comp=[aux_func.dust_vector(data.loc[(data['zobs']>z[0])&(data['zobs']<z[1])].copy(),band_sel,band_1,band_2,
											A_sel,A_b1,A_b2,E_b1b2,mag_cut,b1_cut,b2_cut,
											ebv[(data['zobs']>z[0])&(data['zobs']<z[1])]) for z in zrange]
	
	delta=[]
	mag=[]
	col=[]
	for i in range(len(zrange)):
		delta.append(vector_comp[i]['delta'])
		mag.append(vector_comp[i]['mag'])
		col.append(vector_comp[i]['col'])
		zstr=f'{(zrange[i][0]+zrange[i][1])/2:.2f}'.translate({ord(c): None for c in '.'})
		csv_name=f'{data_dir}dust_vector_{run_name}_z{zstr}_nodust.csv'
		if dusted:
			csv_name=csv_name.replace('nodust','dusted')
		vector_comp[i]['EBV']=EBV
		vector_comp[i].to_csv(csv_name,index=False)
	print('Dust vectors saved')
	
	####################
	# Plots
	####################
	EBV=[EBV]*len(zrange)
	vir=cm.get_cmap('viridis')
	cr=vir(np.array([i/5.0 for i in range(len(zrange))]))
	zlabel=[f'$z_{{{z[0]},{z[1]}}}$' for z in zrange]
	print('Making E(B-V) vs delta/mag/col plot')
	
	plot.figure(figsize=[fs[0],fs[1]*3])
	
	plot.subplot(3,1,1)
	sp.plot(EBV,delta,c=cr,ylabel='$\log(\delta+1)$',title=run_name,plabel=zlabel)
	plot.subplot(3,1,2)
	sp.plot(EBV,mag,c=cr,ylabel=f'${band_1[0]}$ [mag]')
	sp.plot(EBV[0],A_b1,c='k',linestyle='dotted',plabel=f'A$({band_1[0]})$ [mag]')
	plot.subplot(3,1,3)
	sp.plot(EBV,col,c=cr,xlabel='E(B-V) [mag]',ylabel=f'${band_1[0]}-{band_2[0]}$ [mag]')
	sp.plot(EBV[0],E_b1b2,c='k',linestyle='dotted',plabel=f'E$({band_1[0]}-{band_2[0]})$ [mag]')
	plot.tight_layout()
	plot.savefig(f'{plot_dir}dust_vector_{run_name}{fd}.pdf')
	plot.close()
	
	if dusted:
		return(vector_comp)

########################################
# Assign the galaxies to HEALPix pixels
########################################
def pixel_assign(fnames,nside,border_check=False,simple_ebv=True,multithread=False):
	import aux_func
	import numpy as np
	import pandas as pd
	import multiprocessing as mp
	import healpy.pixelfunc as pf
	
	####################
	# Reading data in
	####################
	print('Reading Schlegel map')
	if simple_ebv:
		ebv_map=[np.linspace(0,0.2,pf.nside2npix(n)) for n in nside]
	else:
		ebv_map=[aux_func.sfd_map(percent=False,resample=n,lsst_footprint=False) for n in nside]
	
	print('Reading galaxy data for pixelisation')
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(fnames)))
		temp=[pool.apply_async(aux_func.pix_id,(f,nside,ebv_map,)) for f in fnames]
		temp=[t.get() for t in temp]
	else:
		temp=[aux_func.pix_id(f,nside,ebv_map) for f in fnames]
	res=temp[0][0]
	pix_ids=[t[1] for t in temp]
	pix_ids=[np.concatenate([p[i] for p in pix_ids]) for i in range(len(nside))]
	
	####################
	# Border check
	####################
	print('Checking borders')
	if border_check:
		if multithread:
			pool=mp.Pool(processes=min(multithread,len(fnames)))
			temp=[pool.apply_async(aux_func.find_border,(f,pix_ids,nside,res,)) for f in fnames]
			temp=[t.get() for t in temp]
		else:
			temp=[aux_func.find_border(f,pix_ids,nside,res) for f in fnames]
	
	print('Pixelisation ready')
	return(nside)

########################################
# Calculate pixel properties
########################################
def pixel_stat(fnames,nside,band_sel,band_1,band_2,zrange,mag_cut=24.8,b1_cut=99,b2_cut=99,border_check=False,multithread=False):
	import aux_func
	import numpy as np
	import pandas as pd
	import multiprocessing as mp
	
	####################
	# Stats calculation
	####################
	print('Calculating pixelised properties')
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(fnames)))
		results=[pool.apply_async(aux_func.pix_stat,(f,nside,band_sel,band_1,band_2,mag_cut,b1_cut,b2_cut,zrange,border_check,)) for f in fnames]
		results=[r.get() for r in results]
		pool.close()
	else:
		results=[aux_func.pix_stat(f,nside,band_sel,band_1,band_2,mag_cut,b1_cut,b2_cut,zrange,border_check) for f in fnames]
	results=[r2 for r1 in results for r2 in r1]
	print('Statistics ready')
	return(sorted(results))
