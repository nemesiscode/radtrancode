import numpy as np
from scipy.interpolate import pchip
import matplotlib
matplotlib.use('PDF')
import corner
import matplotlib.pyplot as plt
import pymultinest
import robcat
import nemesisPyMult
import os
import pdb
import pickle
from mpi4py import MPI

# Set plot style
plt.style.use('ggplot')

########################################
# GET MPI INFO (FOR MANY SPECTRAL DRAWS)
########################################

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

###########################
# DEFINE PRIOR + LIKELIHOOD
###########################


def calc_spec(cube, nparams):

# add on H2 and He to VMRs

 vmr = cube[:nvmr]
 vmrsum = np.sum(10**vmr)
 h2 = h2p*(1-vmrsum)
 he = hep*(1-vmrsum)
 vmr = np.append(vmr,[np.log10(h2),np.log10(he)])
 vmr = np.tile(vmr, ( len(pres),1) )

# calculate temperature and height profiles

 mass, rad, T0, alpha, Teff, tau0, n, Pbase, FSH, opac = cube[nvmr:nparams]
 tau0 = 10**tau0
 FSH = 10**FSH
# vmr, mass, rad, T0, alpha, Teff, tau0, n, Pbase, FSH, opac = pickle.load(open('prob.pic','rb'))
# pdb.set_trace()
 (H, temp) = robcat.tempprof(T0,alpha,Teff,tau0,n,vmr,mass,rad)
 if all( H==np.zeros(len(temp)) ):
  return np.zeros(len(spec))

# Find height of cloud

 Hbase= pchip(10**pres[::-1],H[::-1])
 Hb = np.asscalar(Hbase(10**Pbase))

# find folder to calculate in

 ith = 999

# run model
# pickle.dump( (vmr, mass, rad, T0, alpha, Teff, tau0, n, Pbase, FSH, opac) , open('dump.pic', 'wb'))
# pdb.set_trace()

 model = nemesisPyMult.nemesispymult(runname, len(spec), len(temp), vmr.shape[1], ith, temp, vmr, mass, rad, H, Hb, opac, FSH, prad, pvar, pimag, preal, Hb2, opac2, FSH2, nav, flat, flon, solzen, emzen, azi, wt)
# loglikelihood= -0.5*( np.sum( (spec-model)**2/yerr**2 ) )
# pdb.set_trace()
 return model

######################## 
# BEGIN INPUTS
########################

# enter runname
runname = 'gl570d'

# read spectrum and errors
wv, spec, yerr = np.loadtxt("./999/" + runname + ".spx", skiprows=4).T

# define constants, no. gases, h2 percentage, no. params
nvmr=7

h2p = 0.85
hep = 1.0-h2p

n_params=17
nparams=n_params

# give redundant values to unused params

pres=np.log10(np.logspace(3,-4,80))
FSH2=1.0
prad=10.0
pvar=0.1
pimag=1e-5
preal = 1.4
Hb2=10.0
opac2=0.0
nav=0
flat=0.0
flon=0.0
solzen=0.0
emzen=0.0
azi=0.0
wt=0.0


######################## 
# BEGIN FIGURES
########################

# Param names
titles=np.array(['$\log$NH$_3$', '$\log$CH$_4$', '$\log$H$_2$O', '$\log$CO', '$\log$CO$_2$', '$\log$Na','$\log$K', 'Mass [M$_{jup}$]','Radius [R$_{jup}$]', '$T_0$', 'Alpha', '$T_{eff}$', 'Tau_0', '$n$', 'Pbase', 'FSH', 'Opt.Depth'])

'''
#Get retrieval info
prefix='chains/GL570D_'
a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
s=a.get_stats()
values = a.get_equal_weighted_posterior()
'''
values=np.loadtxt('chains/GL570D_phys_live.points')
params=values[:, :n_params]
lnprob=values[:, -1]
samples=params
Nsamp=values.shape[0]
Npars=n_params


# make random draws and plot the spectra
NN=2
draws=np.random.randint(len(samples),size=NN)
xrand=samples[draws,:]
#xrand=samples
model_arr=[]

for i in range(NN):
 model=calc_spec(xrand[i,:], n_params)
 plt.plot(wv,model, lw=1.0, alpha=0.3, color='red')
 model_arr = np.concatenate([model_arr, model])

#plt.semilogy()
plt.errorbar(wv,spec,yerr=yerr,color='black')
plt.xlabel('$\lambda$ ($\mu$m)',size='xx-large')
plt.ylabel('L$_{\lambda}$ [W$\mu$m$^{-1}$]',size='xx-large')
plt.minorticks_on()
plt.tick_params(length=10,width=1,labelsize='xx-large',which='major')
plt.tight_layout()
plt.savefig('Spread_spec.pdf', format='pdf')
plt.close()

# find +/-1 and 2 sigma for spectra
model_arr=model_arr.reshape(NN,wv.shape[0])
pickle.dump(model_arr, open('model_arr.pic','wb') )
mod_median=np.zeros(wv.shape[0])
mod_high_1sig=np.zeros(wv.shape[0])
mod_high_2sig=np.zeros(wv.shape[0])
mod_low_1sig=np.zeros(wv.shape[0])
mod_low_2sig=np.zeros(wv.shape[0])

for i in range(wv.shape[0]):
    percentiles=np.percentile(model_arr[:,i],[4.55, 15.9, 50, 84.1, 95.45])
    mod_low_2sig[i]=percentiles[0]
    mod_low_1sig[i]=percentiles[1]
    mod_median[i]=percentiles[2]
    mod_high_1sig[i]=percentiles[3]
    mod_high_2sig[i]=percentiles[4]

# plot spectral spread
#plt.semilogy()
plt.fill_between(wv,mod_low_2sig,mod_high_2sig,alpha=0.4,edgecolor='None')
plt.fill_between(wv,mod_low_1sig,mod_high_1sig,alpha=0.8,edgecolor='None')
plt.errorbar(wv,spec,yerr=yerr)
plt.plot(wv,mod_median,color='black')

plt.xlabel('$\lambda$ ($\mu$m)',size='xx-large')
plt.ylabel('L$_{\lambda}$ [W$\mu$m$^{-1}$]',size='xx-large')
plt.minorticks_on()
plt.tick_params(length=10,width=1,labelsize='xx-large',which='major')
plt.tight_layout()
plt.savefig('Sum_spec.pdf', format='pdf')
plt.close()

# plot TP profiles

NN=1000
draws=np.random.randint(len(samples),size=NN)
xrand=samples[draws,:]
Tarr=[]
P = np.logspace(np.log10(1000.0),-4.0,80)
for i in range(NN):
 cube = xrand[i,:]
 vmr = cube[:nvmr]
 vmrsum = np.sum(10**vmr)
 h2 = h2p*(1-vmrsum)
 he = hep*(1-vmrsum)
 vmr = np.append(vmr,[np.log10(h2),np.log10(he)])
 vmr = np.tile(vmr, ( len(pres),1) )
 mass, rad, T0, alpha, Teff, tau0, n, Pbase, FSH, opac = cube[nvmr:nparams]
 tau0 = 10**tau0
 (H, temp) = robcat.tempprof(T0,alpha,Teff,tau0,n,vmr,mass,rad)
 Tarr=np.concatenate([Tarr,temp])

Tarr=Tarr.reshape(NN,P.shape[0])

for i in range(NN):
 plt.plot(Tarr[i,:],P, lw=0.5, alpha=0.1, color='red')

plt.minorticks_on()
plt.tick_params(length=10,width=1,labelsize='xx-large',which='major')
plt.xlabel('Temperature [K]',size='xx-large')
plt.ylabel('Pressure [bar]',size='xx-large')
plt.semilogy()
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('Spag_TP.pdf',format='pdf')
plt.close()

#plotting spread
Tmedian=np.zeros(P.shape[0])
Tlow_1sig=np.zeros(P.shape[0])
Thigh_1sig=np.zeros(P.shape[0])
Tlow_2sig=np.zeros(P.shape[0])
Thigh_2sig=np.zeros(P.shape[0])

for i in range(P.shape[0]):
    percentiles=np.percentile(Tarr[:,i],[4.55, 15.9, 50, 84.1, 95.45])
    Tlow_2sig[i]=percentiles[0]
    Tlow_1sig[i]=percentiles[1]
    Tmedian[i]=percentiles[2]
    Thigh_1sig[i]=percentiles[3]
    Thigh_2sig[i]=percentiles[4]

plt.fill_betweenx(P,Tlow_2sig,Thigh_2sig,facecolor='r',edgecolor='None',alpha=0.3)
plt.fill_betweenx(P,Tlow_1sig,Thigh_1sig,facecolor='r',edgecolor='None',alpha=1.)
plt.plot(Tmedian, P,'b')
plt.semilogy()
plt.gca().invert_yaxis()
plt.xlabel('Temperature [K]',size='xx-large')
plt.ylabel('Pressure [bar]',size='xx-large')
plt.minorticks_on()
plt.tick_params(length=10,width=1,labelsize='xx-large',which='major')
plt.tight_layout()
plt.savefig('Sum_TP.pdf',format='pdf')
plt.close()

# plot samples

#matplotlib.rcParams.update({'font.size': 18})
figure=corner.corner(samples, bin=200, quantiles=[0.16,0.5,0.84], labels=titles, show_titles='True', plot_contours='True')
figure.savefig("triangle.pdf")
