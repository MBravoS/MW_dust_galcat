import glob
import numpy as np
import pandas as pd
import splotch as sp
import cmocean as cmo
import multiprocessing as mp
import matplotlib.pyplot as plot
import matplotlib.gridspec as gs

sp.use_style('/home/mbravo/pypati.style')
new_fig_size=np.array(plot.rcParams.get('figure.figsize')) 
########################################
# Loading data
########################################
def csv_load(f,m1,m2,m3):
	df=pd.read_csv(f)
	df=df.loc[(df['r_ap_nodust']<m1)&(df['u_ap_nodust']<m2)&(df['z_ap_nodust']<m3)]
	return(df)

GAL248_files=glob.glob('/fast_scratch2/mbravo/MWdust_data/galaxies*GALFORM*csv')
Buz248_files=glob.glob('/fast_scratch2/mbravo/MWdust_data/galaxies*Buzzard*csv')
Buz260_files=glob.glob('/fast_scratch2/mbravo/MWdust_data_deep/galaxies*Buzzard*csv')

pool=mp.Pool(processes=min(56,len(GAL248_files)))
temp=[pool.apply_async(csv_load,(f,24.8,27.3,27.2,)) for f in GAL248_files]
temp=[t.get() for t in temp]
GAL248_data=pd.concat(temp)
pool=mp.Pool(processes=min(56,len(Buz248_files)))
temp=[pool.apply_async(csv_load,(f,24.8,27.3,27.2,)) for f in Buz248_files]
temp=[t.get() for t in temp]
Buz248_data=pd.concat(temp)
pool=mp.Pool(processes=min(56,len(Buz260_files)))
temp=[pool.apply_async(csv_load,(f,24.8,27.3,27.2,)) for f in Buz260_files]
temp=[t.get() for t in temp]
Buz260_data=pd.concat(temp)

zrange=[[0.0,0.3],[0.3,0.6],[0.6,0.9],[0.9,1.2],[1.2,2.5]]

########################################
# Plots
########################################
colour_list=[]
for i in np.linspace(0,1,5):
	temp_col=[np.array(cmo.cm.haline(i)) for j in range(2)]
	temp_col[0][-1]=0.4
	temp_col[1][-1]=0.8
	colour_list.append(temp_col)

####################
# col-mag (ap)
####################
fig=plot.figure(figsize=(new_fig_size[0],new_fig_size[0]*2))
spec=gs.GridSpec(nrows=3,ncols=1,figure=fig,wspace=0,hspace=0,left=0.17,right=0.98,bottom=0.08,top=0.97)
fax=[fig.add_subplot(spec[0,0]),fig.add_subplot(spec[1,0]),fig.add_subplot(spec[2,0])]
temp=fax[0].get_xaxis().set_ticklabels([])
temp=fax[1].get_xaxis().set_ticklabels([])
for i in range(3):
	#temp=fax[i].get_xaxis().set_ticks([0,5e3,10e3,15e3,20e3,25e3])
	#temp=fax[i].get_yaxis().set_ticks([0,15e-5,30e-5,45e-5,60e-5,75e-5,90e-5])
	#temp=fax[i].get_yaxis().set_ticklabels(['0','$1.5\\textsc{e}^{-3}$','$3.0\\textsc{e}^{-3}$','$4.5\\textsc{e}^{-3}$',
	#											'$6.0\\textsc{e}^{-3}$','$7.5\\textsc{e}^{-3}$','$9.0\\textsc{e}^{-3}$'])
	fax[i].set_rasterized(True)
j=0 
for zr in zrange:
	if j%2==0:
		temp=GAL248_data.loc[(GAL248_data['zobs_sim']>zr[0])&(GAL248_data['zobs_sim']<zr[1])]
		sp.contourp(temp['u_ap_nodust'],temp['u_ap_nodust']-temp['z_ap_nodust'],bins=[np.linspace(19,28,61),np.linspace(-1,7,31)],ax=fax[0],
					xinvert=True,smooth=0.8,filled=True,colors=colour_list[j],plabel=False,xlabel='$u$ [mag]',ylabel='$u-z$ [mag]')
		temp=Buz248_data.loc[(Buz248_data['zobs_sim']>zr[0])&(Buz248_data['zobs_sim']<zr[1])]
		sp.contourp(temp['u_ap_nodust'],temp['u_ap_nodust']-temp['z_ap_nodust'],bins=[np.linspace(19,28,61),np.linspace(-1,7,31)],ax=fax[1],
					xinvert=True,smooth=0.8,filled=True,colors=colour_list[j],plabel=False,xlabel='$u$ [mag]',ylabel='$u-z$ [mag]')
		temp=Buz260_data.loc[(Buz260_data['zobs_sim']>zr[0])&(Buz260_data['zobs_sim']<zr[1])]
		sp.contourp(temp['u_ap_nodust'],temp['u_ap_nodust']-temp['z_ap_nodust'],bins=[np.linspace(19,28,61),np.linspace(-1,7,31)],ax=fax[2],
					xinvert=True,smooth=0.8,filled=True,colors=colour_list[j],plabel=False,xlabel='$u$ [mag]',ylabel='$u-z$ [mag]')
	j+=1
plot.savefig('/fast_scratch2/mbravo/MWdust_plots/mag_col_ap.pdf')
plot.savefig('/fast_scratch2/mbravo/MWdust_plots/mag_col_ap.png')
plot.close()

####################
# col-mag (ab)
####################
fig=plot.figure(figsize=(new_fig_size[0],new_fig_size[0]*2))
spec=gs.GridSpec(nrows=3,ncols=1,figure=fig,wspace=0,hspace=0,left=0.17,right=0.98,bottom=0.08,top=0.97)
fax=[fig.add_subplot(spec[0,0]),fig.add_subplot(spec[1,0]),fig.add_subplot(spec[2,0])]
temp=fax[0].get_xaxis().set_ticklabels([])
temp=fax[1].get_xaxis().set_ticklabels([])
for i in range(3):
	#temp=fax[i].get_xaxis().set_ticks([0,5e3,10e3,15e3,20e3,25e3])
	#temp=fax[i].get_yaxis().set_ticks([0,15e-5,30e-5,45e-5,60e-5,75e-5,90e-5])
	#temp=fax[i].get_yaxis().set_ticklabels(['0','$1.5\\textsc{e}^{-3}$','$3.0\\textsc{e}^{-3}$','$4.5\\textsc{e}^{-3}$',
	# 											'$6.0\\textsc{e}^{-3}$','$7.5\\textsc{e}^{-3}$','$9.0\\textsc{e}^{-3}$'])
	fax[i].set_rasterized(True) 
j=0
for zr in zrange:
	if j%2==0:
		temp=GAL248_data.loc[(GAL248_data['zobs_sim']>zr[0])&(GAL248_data['zobs_sim']<zr[1])]
		sp.contourp(temp['u_ab'],temp['u_ab']-temp['z_ab'],bins=[np.linspace(-22.5,-10.5,31),np.linspace(-1,4,31)],
					ax=fax[0],xinvert=True,smooth=0.8,filled=True,colors=colour_list[j],plabel=False,
					ylabel='${u-z}_\mathrm{ab}$ [mag]')
		temp=Buz248_data.loc[(Buz248_data['zobs_sim']>zr[0])&(Buz248_data['zobs_sim']<zr[1])]
		sp.contourp(temp['u_ab'],temp['u_ab']-temp['z_ab'],bins=[np.linspace(-22.5,-10.5,31),np.linspace(-1,4,31)],
					ax=fax[1],xinvert=True,smooth=0.8,filled=True,colors=colour_list[j],plabel=False,
					ylabel='${u-z}_\mathrm{ab}$ [mag]')
		temp=deep_Buz260_data.loc[(deep_Buz260_data['zobs_sim']>zr[0])&(deep_Buz260_data['zobs_sim']<zr[1])]
		sp.contourp(temp['u_ab'],temp['u_ab']-temp['z_ab'],bins=[np.linspace(-22.5,-10.5,31),np.linspace(-1,4,31)],
					ax=fax[2],xinvert=True,smooth=0.8,filled=True,colors=colour_list[j],plabel=False,
					xlabel='$u_\mathrm{ab}$ [mag]',ylabel='${u-z}_\mathrm{ab}$ [mag]')
	j+=1
plot.savefig('/fast_scratch2/mbravo/MWdust_plots/mag_col_ab.pdf')
plot.savefig('/fast_scratch2/mbravo/MWdust_plots/mag_col_ab.png')
plot.close()

####################
# col-col
####################
fig=plot.figure(figsize=(new_fig_size[0],new_fig_size[0]*2))
spec=gs.GridSpec(nrows=3,ncols=1,figure=fig,wspace=0,hspace=0,left=0.17,right=0.98,bottom=0.08,top=0.97)
fax=[fig.add_subplot(spec[0,0]),fig.add_subplot(spec[1,0]),fig.add_subplot(spec[2,0])]
temp=fax[0].get_xaxis().set_ticklabels([])
temp=fax[1].get_xaxis().set_ticklabels([])
for i in range(3):
	#temp=fax[i].get_xaxis().set_ticks([0,5e3,10e3,15e3,20e3,25e3])
	#temp=fax[i].get_yaxis().set_ticks([0,15e-5,30e-5,45e-5,60e-5,75e-5,90e-5])
	#temp=fax[i].get_yaxis().set_ticklabels(['0','$1.5\\textsc{e}^{-3}$','$3.0\\textsc{e}^{-3}$','$4.5\\textsc{e}^{-3}$',
	#											'$6.0\\textsc{e}^{-3}$','$7.5\\textsc{e}^{-3}$','$9.0\\textsc{e}^{-3}$'])
	fax[i].set_rasterized(True)
temp=Buz260_data.loc[(Buz260_data['zobs_sim']>zrange[0][0])&(Buz260_data['zobs_sim']<zrange[0][1])]
sp.hist2D(temp['r_ab']-temp['z_ab'],temp['u_ab']-temp['r_ab'],bins=[np.linspace(-1,1.5,31),np.linspace(-1,3,31)],
			ax=fax[0],cmap=cmo.cm.gray_r,clog=False,ylabel='${u-r}_\mathrm{ab}$ [mag]')
temp=Buz260_data.loc[(Buz260_data['zobs_sim']>zrange[2][0])&(Buz260_data['zobs_sim']<zrange[2][1])]
sp.hist2D(temp['r_ab']-temp['z_ab'],temp['u_ab']-temp['r_ab'],bins=[np.linspace(-1,1.5,31),np.linspace(-1,3,31)],
			ax=fax[1],cmap=cmo.cm.gray_r,ylabel='${u-r}_\mathrm{ab}$ [mag]')
temp=Buz260_data.loc[(Buz260_data['zobs_sim']>zrange[4][0])&(Buz260_data['zobs_sim']<zrange[4][1])]
sp.hist2D(temp['r_ab']-temp['z_ab'],temp['u_ab']-temp['r_ab'],bins=[np.linspace(-1,1.5,31),np.linspace(-1,3,31)],
			ax=fax[2],cmap=cmo.cm.gray_r,xlabel='${r-z}_\mathrm{ab}$ [mag]',ylabel='${u-r}_\mathrm{ab}$ [mag]')
plot.savefig('/fast_scratch2/mbravo/MWdust_plots/col_col.pdf')
plot.savefig('/fast_scratch2/mbravo/MWdust_plots/col_col.png')
plot.close()
