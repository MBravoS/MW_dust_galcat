'''Here are the core functions for the map analysis'''
########################################
# Dust map recovery
########################################
def dust_mapping(pnames,dvec,nside,zrange,out_dir,plot_dir):
	import time
	import aux_func
	import numpy as np
	import pandas as pd
	import splotch as sp
	import scipy.optimize as opti
	import matplotlib.lines as lines
	import matplotlib.pyplot as plot
	import scipy.interpolate as interp
	import matplotlib.patches as patches
	
	t0=time.time()
	for i in range(len(nside)):
		ns=nside[i]
		nside_key=f'n{np.log2(ns):.0f}'
		data=pd.concat([pd.read_csv(p) for p in pnames if nside_key in p],ignore_index=True)
		Debv,ebv_from_delta,ebv_from_mag,ebv_from_col=[],[],[],[]
		EBV_recovery,z_label=[],[]
		EBV_input=data[f'{nside_key}_EBV']
		for j in range(len(zrange)):
			zr=zrange[j]
			z_label.append(f'$z_{{{np.sum(zr)/2:.2f}}}$')
			z_key=f'z{np.sum(zr)/2:.2f}'.replace('.','')
			
			####################
			# Delta calc
			####################
			delta=data[f'{nside_key}_{z_key}_count']+np.random.uniform(-0.5,0.5,len(data))
			data[f'{nside_key}_{z_key}_delta']=np.log10(delta/np.mean(delta))
			
			####################
			# Read dust vector
			####################
			dust_vector=pd.read_csv(f'{out_dir}dust_vector_{pnames[0].split("/")[-1].split("_")[2]}_{z_key}.csv')
			dust_vector,ebv2d,ebv2m,ebv2c=aux_func.slope2(dust_vector,np.array(data[f'{nside_key}_EBV']))
			
			plot.figure()
			sp.plot(dust_vector['deltaEBV'],dust_vector['delta'],plabel='$\log(\delta+1)_\mathrm{measured}$')
			sp.plot(dust_vector['deltaEBV'],dust_vector['mag'],plabel='$u_\mathrm{measured}$')
			sp.plot(dust_vector['deltaEBV'],dust_vector['col'],plabel='$(u-z)_\mathrm{measured}$')
			sp.plot(dust_vector['deltaEBV'],dust_vector['deltaEBV']*ebv2d,linestyle='dotted',plabel='$\delta_\mathrm{fit}$')
			sp.plot(dust_vector['deltaEBV'],dust_vector['deltaEBV']*ebv2m,linestyle='dotted',plabel='$u_\mathrm{fit}$')
			sp.plot(dust_vector['deltaEBV'],dust_vector['deltaEBV']*ebv2c,linestyle='dotted',plabel='${u-z}_\mathrm{fit}$',
					xlabel='$\Delta E(B-V)$ [mag]',ylabel='$\Delta$',xlim=[-0.1,0.1])
			plot.tight_layout()
			plot.savefig(f'{plot_dir}metric_check_{pnames[0].split("/")[-1].split("_")[2]}_{z_key}.pdf')
			plot.close()
			
			Debv.append(data[f'{nside_key}_EBV']-np.median(data[f'{nside_key}_EBV']))
			#ebv_from_delta.append(data[f'{nside_key}_{z_key}_delta']/ebv2d)
			ebv_from_delta.append((data[f'{nside_key}_{z_key}_delta']-np.median(data[f'{nside_key}_{z_key}_delta']))/ebv2d)
			ebv_from_mag.append((data[f'{nside_key}_{z_key}_mag']-np.median(data[f'{nside_key}_{z_key}_mag']))/ebv2m)
			ebv_from_col.append((data[f'{nside_key}_{z_key}_col']-np.median(data[f'{nside_key}_{z_key}_col']))/ebv2c)
			
			spl_delta=interp.UnivariateSpline(dust_vector['deltaEBV'],dust_vector['delta'],ext=0)
			spl_mag=interp.UnivariateSpline(dust_vector['deltaEBV'],dust_vector['mag'],ext=0)
			spl_col=interp.UnivariateSpline(dust_vector['deltaEBV'],dust_vector['col'],ext=0)
			
			def dist(val,D,M,C):
				dd=spl_delta(val)
				mm=spl_mag(val)
				cc=spl_col(val)
				return(((D-dd/ebv2d)**2+(M-mm/ebv2m)**2+(C-cc/ebv2c)**2)**0.5)
			
			ebv_recover=np.zeros(len(ebv_from_delta[j]))
			for k in range(len(ebv_from_delta[j])):
				ebv_recover[k]=opti.minimize(dist,x0=0.0,args=(ebv_from_delta[j][k],ebv_from_mag[j][k],ebv_from_col[j][k])).x
			
			EBV_recovery.append(ebv_recover+np.median(EBV_input))
		
		EBV_final=np.array(EBV_recovery)
		EBV_final=np.average(EBV_final,axis=0,weights=1/np.var(EBV_final,axis=1))
		
		####################
		# Plots
		####################
		msize=20
		if nside[i]>65:
			msize=1
		
		plot.figure()
		for j in range(len(zrange)):
			sp.scatter(EBV_input-np.median(EBV_input),ebv_from_delta[j],c=f'C{j}',marker='d',s=msize)
			sp.scatter(EBV_input-np.median(EBV_input),ebv_from_mag[j],c=f'C{j}',marker='*',s=msize)
			sp.scatter(EBV_input-np.median(EBV_input),ebv_from_col[j],c=f'C{j}',xlabel='$\Delta E(B-V)_\mathrm{input}$',
						ylabel='$\Delta E(B-V)_\mathrm{recovered}$',xlim=[-0.1,0.2],s=msize)
		sp.axline(a=1,color='k',linestyle='dashed')
		plot.legend([patches.Patch(color=f'C{j}') for j in range(len(zrange))]+
					[lines.Line2D([0.5],[0.5],color='gray',marker='d',linestyle=''),
						lines.Line2D([0.5],[0.5],color='gray',marker='*',linestyle=''),
						lines.Line2D([0.5],[0.5],color='gray',marker='o',linestyle=''),
						lines.Line2D([0,1],[0,1],color='black',linestyle='dashed')],
					[f'{z_label[j]}' for j in range(len(zrange))]+['$\delta$','$u$','$u-z$','1:1'])
		plot.tight_layout()
		plot.savefig(f'{plot_dir}delta_ebv_recovery_zbin_{pnames[0].split("/")[-1].split("_")[2]}_{nside_key}_full.pdf')
		
		plot.figure()
		for j in range(len(zrange)):
			sp.scatter(EBV_input,EBV_recovery[j],plabel=z_label[j],c=f'C{j}',xlabel='$E(B-V)_\mathrm{input}$',
						ylabel='$E(B-V)_\mathrm{recovered}$',xlim=[0,0.2],ylim=[-0.1,0.3],s=msize)
		sp.axline(a=1,plabel='1:1',color='k',linestyle='dashed')
		plot.tight_layout()
		plot.savefig(f'{plot_dir}ebv_recovery_zbin_{pnames[0].split("/")[-1].split("_")[2]}_{nside_key}_combined.pdf')
		
		plot.figure()
		sp.scatter(EBV_input,EBV_final,c='C0',xlabel='$E(B-V)_\mathrm{input}$',ylabel='$E(B-V)_\mathrm{recovered}$',
					xlim=[0,0.2],ylim=[0,0.2],s=msize)
		sp.axline(a=1,plabel='1:1',color='k',linestyle='dashed')
		plot.tight_layout()
		plot.savefig(f'{plot_dir}ebv_recovery_zbin_{pnames[0].split("/")[-1].split("_")[2]}_{nside_key}_final.pdf')
	
	t1=time.time()
	print(f'Dust map recovered in {t1-t0} s')

#def dust_mapping2(pnames,dvec,nside,zrange,out_dir,plot_dir):
#	import aux_func
#	import numpy as np
#	import pandas as pd
#	import splotch as sp
#	import scipy.optimize as opti
#	import matplotlib.pyplot as plot
#	import scipy.interpolate as interp
#	
#	for i in range(len(nside)):
#		ns=nside[i]
#		nside_key=f'n{np.log2(ns):.0f}'
#		data=pd.concat([pd.read_csv(p) for p in pnames if nside_key in p],ignore_index=True)
#		EBV_recovery=[]
#		z_label=[]
#		for j in range(len(zrange)):
#			zr=zrange[j]
#			z_key=f'z{np.sum(zr)/2:.2f}'.replace('.','')
#			
#			####################
#			# Delta calc
#			####################
#			delta=data[f'{nside_key}_{z_key}_count']+np.random.uniform(-0.5,0.5,len(data))
#			data[f'{nside_key}_{z_key}_delta']=np.log10(delta/np.mean(delta))
#			
#			####################
#			# Read dust vector
#			####################
#			dust_vector=pd.read_csv(f'{out_dir}dust_vector_{pnames[0].split("/")[-1].split("_")[2]}_{z_key}_dusted.csv')
#			dust_vector,ebv2d,ebv2m,ebv2c=aux_func.slope2(dust_vector,np.array(data[f'{nside_key}_EBV']))
#			
#			ebv_from_delta=(data[f'{nside_key}_{z_key}_delta']-np.median(data[f'{nside_key}_{z_key}_delta']))/ebv2d
#			ebv_from_mag=(data[f'{nside_key}_{z_key}_mag']-np.median(data[f'{nside_key}_{z_key}_mag']))/ebv2m
#			ebv_from_col=(data[f'{nside_key}_{z_key}_col']-np.median(data[f'{nside_key}_{z_key}_col']))/ebv2c
#			
#			spl_delta=interp.UnivariateSpline(dust_vector['deltaEBV'],dust_vector['delta'])
#			spl_mag=interp.UnivariateSpline(dust_vector['deltaEBV'],dust_vector['mag'])
#			spl_col=interp.UnivariateSpline(dust_vector['deltaEBV'],dust_vector['col'])
#			
#			#sp.plot([dust_vector['deltaEBV']]*3,
#			#		[spl_delta(dust_vector['deltaEBV']),spl_mag(dust_vector['deltaEBV']),spl_col(dust_vector['deltaEBV'])],
#			#		plabel=['d','m','c'])
#			#plot.show()
#			
#			def dist(val,D,M,C):
#				dd=spl_delta(val)
#				mm=spl_mag(val)
#				cc=spl_col(val)
#				return(((D-dd)**2+(M-mm)**2+(C-cc)**2)**0.5)
#			
#			ebv_recover=np.zeros(len(ebv_from_delta))
#			for k in range(len(ebv_from_delta)):
#				ebv_recover[k]=opti.minimize(dist,x0=0.0,args=(ebv_from_delta[k],ebv_from_mag[k],ebv_from_col[k])).x
#			
#			EBV_recovery.append(ebv_recover)
#			z_label.append(z_key)
#		
#		####################
#		# Plots
#		####################
#		plot.figure()
#		
#		for j in range(len(zrange)):
#			sp.scatter(data[f'{nside_key}_EBV']-np.median(data[f'{nside_key}_EBV']),EBV_recovery[j]+0.1,plabel=z_label[j],c=f'C{j}',
#						xlabel='$\Delta E(B-V)_\mathrm{input}$',ylabel='$\Delta E(B-V)_\mathrm{recovered}$')
#		sp.axline(m=1,plabel='1:1',color='k',linestyle='dashed')
#		plot.tight_layout()
#		plot.savefig(f'{plot_dir}delta_ebv_recovery_zbin_{pnames[0].split("/")[-1].split("_")[2]}.pdf')

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
	fd=''
	
	####################
	# Reading data in
	####################
	t0=time.time()
	print('Reading galaxy data for dust vector calculation')
	run_name=fnames[0].split("galaxies_")[1].split('_')[0]
	data=[]
	for f in fnames:
		print(f)
		data.append(pd.read_csv(f))
	data=pd.concat(data)
	ebv=np.zeros(len(data))
	
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
		vector_comp=[pool.apply(aux_func.dust_vector,(data.loc[(data['zobs']>z[0])&(data['zobs']<z[1])].copy(),band_sel,band_1,band_2,
														EBV,mag_cut,b1_cut,b2_cut, ebv[(data['zobs']>z[0])&(data['zobs']<z[1])],)) for z in zrange]
	else:
		vector_comp=[aux_func.dust_vector(data.loc[(data['zobs']>z[0])&(data['zobs']<z[1])].copy(),band_sel,band_1,band_2,
											EBV,mag_cut,b1_cut,b2_cut,ebv[(data['zobs']>z[0])&(data['zobs']<z[1])]) for z in zrange]
	#vector_comp=[aux_func.dust_vector(data.loc[(data['zobs']>z[0])&(data['zobs']<z[1])].copy(),band_sel,band_1,band_2,
	#									EBV,mag_cut,b1_cut,b2_cut,ebv[(data['zobs']>z[0])&(data['zobs']<z[1])]) for z in zrange]
	
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
	EBV=[EBV]*len(zrange)
	vir=cm.get_cmap('viridis')
	cr=vir(np.array([i/5.0 for i in range(len(zrange))]))
	zlabel=[f'$z_{{{z[0]},{z[1]}}}$' for z in zrange]
	
	print('Making E(B-V) vs delta/mag/col plot')
	
	plot.figure(figsize=[fs[0],fs[1]*3])
	
	plot.subplot(3,1,1)
	sp.plot(EBV,delta,c=cr,ylabel='$\Delta\log(\delta+1)$',title=run_name,plabel=zlabel)
	plot.subplot(3,1,2)
	sp.plot(EBV,mag,c=cr,ylabel=f'$\Delta {band_1[0]}$ [mag]')
	sp.plot(EBV[0],A_b1,c='k',linestyle='dashed',plabel=f'A$({band_1[0]})$ [mag]')
	plot.subplot(3,1,3)
	sp.plot(EBV,col,c=cr,xlabel='E(B-V) [mag]',ylabel=f'$\Delta({band_1[0]}-{band_2[0]})$ [mag]')
	sp.plot(EBV[0],E_b1b2,c='k',linestyle='dashed',plabel=f'E$({band_1[0]}-{band_2[0]})$ [mag]')
	plot.tight_layout()
	plot.savefig(f'{plot_dir}dust_vector_{run_name}{fd}.pdf')
	plot.close()
	
	t1=time.time()
	print(f'Dust vectors created in {t1-t0} s')
	
	return(vector_comp)

########################################
# Add observational errors to galaxies
########################################
def magz_err(fnames,multithread):
	import aux_func
	import multiprocessing as mp
	
	if multithread:
		pool=mp.Pool(processes=min(multithread,len(fnames)))
		temp=[pool.apply_async(aux_func.magz_err_perfile,(f,)) for f in fnames]
		temp=[t.get() for t in temp]
	else:
		temp=[aux_func.magz_err_perfile(f) for f in fnames]

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
	#temp=[aux_func.pix_id(fnames[j],nside,ebv_map,j) for j in range(len(fnames))]
	
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
def pixel_stat(fnames,nside,band_sel,band_1,band_2,zrange,mag_cut=24.8,b1_cut=99,b2_cut=99,border_check=False,multithread=False):
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
		results=[pool.apply_async(aux_func.pix_stat,(f,nside,band_sel,band_1,band_2,mag_cut,b1_cut,b2_cut,zrange,border_check,)) for f in fnames]
		results=[r.get() for r in results]
		pool.close()
	else:
		results=[aux_func.pix_stat(f,nside,band_sel,band_1,band_2,mag_cut,b1_cut,b2_cut,zrange,border_check) for f in fnames]
	#results=[aux_func.pix_stat(f,nside,band_sel,band_1,band_2,mag_cut,b1_cut,b2_cut,zrange,border_check) for f in fnames]
	
	results=[r2 for r1 in results for r2 in r1 if r2 is not None]
	
	t1=time.time()
	print(f'Pixel statistics finished in {t1-t0} s')
	return(sorted(results))
