'''Here are the core functions for the map analysis'''

########################################
# Data characterization
########################################
def data_char(fnames,band_sel,band_1,band_2,data_dir,plot_dir,z_range=[[0.0,0.3],[0.3,0.6],[0.6,0.9],[0.9,1.2],[1.2,2.5]],
				mag_cut=24.8,sigma_cut=False,multithread=False):
	import aux_func
	import numpy as np
	import pandas as pd
	import splotch as sp
	#import astropy.units as u
	import matplotlib.cm as cm
	import multiprocessing as mp
	import matplotlib.pyplot as plot
	#from astropy.cosmology import WMAP7 as cosmo
	sp.use_style('splotch.style')
	fs=np.array(plot.rcParams.get('figure.figsize'))
	
	####################
	# Reading data in
	####################
	data=pd.concat([pd.read_csv(f) for f in fnames])
	
	####################
	# Extinction effect
	####################
	EBV=np.linspace(0,0.2,num=200)
	A_sel,temp=aux_func.extinction_law(EBV,band_sel,band_1)
	A_b1,E_b1b2=aux_func.extinction_law(EBV,band_1,band_2)
	A_b2=A_b1-E_b1b2
	
	if multithread:
		cpus=min(multithread,len(Ar))
		pool1=mp.Pool(processes=cpus)
		vector_comp=[pool1.apply(aux_func.dust_vector(data,band_sel,band_1,band_2,z,A_sel,A_b1,E_b1b2,mag_cut,)) for z in z_range]
	else:
		vector_comp=[aux_func.dust_vector(data,band_sel,band_1,band_2,z,A_sel,A_b1,E_b1b2,mag_cut) for z in z_range]
	
	delta=[]
	mag=[]
	col=[]
	for i in range(len(z_range)):
		delta.append(vector_comp[i]['delta'])
		mag.append(vector_comp[i]['mag'])
		col.append(vector_comp[i]['col'])
		vector_comp[i].to_csv(f'{data_dir}dust_vector_comp_z_{(z_range[i][0]+z_range[i][1])/2:.1f}.csv',index=False)
	
	####################
	# Plots
	####################
	EBV=[EBV]*len(z_range)
	vir=cm.get_cmap('viridis')
	cr=vir(np.array([i/5.0 for i in range(len(z_range))]))
	zlabel=[f'$z_{{{z[0]},{z[1]}}}$' for z in z_range]
	
	plot.figure(figsize=[fs[0],fs[1]*3])
	
	plot.subplot(3,1,1)
	sp.plot(EBV,delta,c=cr,ylabel='$\log(\delta+1)$',title=fnames[0].split('/')[-1].split('_galaxies_')[0],plabel=zlabel)
	plot.subplot(3,1,2)
	sp.plot(EBV,mag,c=cr,ylabel=f'${band_1[0]}$ [mag]')
	sp.plot(EBV[0],A_b1,c='k',linestyle='dotted',plabel=f'A$({band_1[0]})$ [mag]')
	plot.subplot(3,1,3)
	sp.plot(EBV,col,c=cr,xlabel='E(B-V) [mag]',ylabel=f'${band_1[0]}-{band_2[0]}$ [mag]')
	sp.plot(EBV[0],E_b1b2,c='k',linestyle='dotted',plabel=f'E$({band_1[0]}-{band_2[0]})$ [mag]')
	plot.tight_layout()
	plot.savefig(f'{plot_dir}dust_vector.pdf')
	plot.close()

########################################
# Assign the galaxies to HEALPix pixels
########################################
def pixel_assign(fnames,nside=[64],border_check=False,multithread=False):
	import aux_func
	import numpy as np
	import pandas as pd
	import healpy as hp
	import multiprocessing as mp
	
	####################
	# Reading data in
	####################
	if multithread:
		cpus=min(multithread,len(fnames))
		pool1=mp.Pool(processes=cpus)
		temp=[pool1.apply_async(aux_func.pix_id,(f,nside,)) for f in fnames]
		temp=[t.get() for t in temp]
	else:
		temp=[aux_func.pix_id(f,nside) for f in fnames]
	res=temp[0][0]
	pix_ids=[t[1] for t in temp]
	pix_ids=[np.concatenate([p[i] for p in pix_ids]) for i in range(len(nside))]
	
	####################
	# Border check
	####################
	if border_check:
		if multithread:
			cpus=min(multithread,len(fnames))
			pool2=mp.Pool(processes=cpus)
			temp=[pool2.apply_async(aux_func.find_border,(f,pix_ids,nside,res,)) for f in fnames]
			temp=[t.get() for t in temp]
		else:
			temp=[aux_func.find_border(f,pix_ids,nside,res) for f in fnames]
	
	return(nside)

########################################
# Calculate pixel properties
########################################
#def pixel_stats(fnames,nside,band_1,band_2,cc_bands,redshift_cuts=[0.3,0.6,0.9,1.2],mag_cut=False,sigma_cut=False):
#	import numpy as np
#	import pandas as pd
#	import multiprocessing as mp
#	
#	Lz=len(redshift_cuts)+1
#	Ln=len(nside)
#	########################
#	#Stats calculation
#	
#	EBV_map=[base.Map_Analysis(percent=False,resample=n) for n in nside]
#	pool2=mp.Pool(processes=cpus)
#	results2=[[[pool2.apply_async(base.Pix_Stats_Calc,(nside[j],fnames[j][k],zcond[i][j][k],band_1,band_2,cc_bands,r_band,
#									r_cut,EBV_map[j],sigma_cut,)) for k in xrange(len(fnames[j]))] for j in xrange(L2)] for i in xrange(L1)]
#	results2=[[[p3.get() for p3 in p2] for p2 in p1] for p1 in results2]
#	pool2.close()
#	#results2=[[[base.Pix_Stats_Calc(nside[j],fnames[j][k],zcond[i][j][k],band_1,band_2,cc_bands,r_band,
#	#								r_cut,EBV_map[j],sigma_cut) for k in xrange(len(fnames[j]))] for j in xrange(L2)] for i in xrange(L1)]
#	print 'Basic stats ready'
#	stats_dict={}
#	#Pix values
#	stats_dict['counts']=[[np.concatenate([fn[0] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['mag']=[[np.concatenate([fn[1] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['colour']=[[np.concatenate([fn[2] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['redshift']=[[np.concatenate([fn[3] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['colour_comp_1']=[[np.concatenate([fn[4] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['colour_comp_2']=[[np.concatenate([fn[5] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['pix_ID']=[[np.concatenate([fn[6] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['counts_AE']=[[np.concatenate([fn[7] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['mag_A']=[[np.concatenate([fn[8] for fn in ns]) for ns in zc] for zc in results2]
#	stats_dict['col_E']=[[np.concatenate([fn[9] for fn in ns]) for ns in zc] for zc in results2]
#	del results2
#	with open(data_dir+lc+mag_type+'_stats.dat', 'wb') as fp:
#		pickle.dump(stats_dict,fp)
