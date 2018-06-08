import numpy as np
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
from scipy.interpolate import pchip
from pylab import plot, show


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

def t_curve(A, B, P):
 tcon = 1e4 / (A + B*np.log10(P) )
 return tcon

def svp_curve(A,B,T):
 svp = 10**(A + B/T)
 return svp

def calc_spec(cube, nparams, ith):

# add on H2 and He to VMRs
 P = np.logspace(3,-2,80)
 vmr = cube[:nvmr]
 vmrsum = np.sum(10**vmr)
 h2 = h2p*(1-vmrsum)
 he = hep*(1-vmrsum)
 vmr = np.append(vmr,[np.log10(h2),np.log10(he)])
 vmr = np.tile(vmr, ( len(P),1) )

# calculate temperature and height profiles

 mass, rad, T0, alpha, Teff, tau0, n, Pbase, FSH, opac = cube[nvmr:nparams]
 FSH= 10**FSH
 tau0 = 10**tau0
 pc = 3.086e+16

# pickle.dump((vmr,mass, rad, T0, alpha, Teff, tau0, n), open(str(ith)+'.pic','wb'))

 (H, temp) = robcat.tempprof(T0,alpha,Teff,tau0,n,vmr,mass,rad)
 if all(H==np.zeros(len(temp)) ):
  return -np.inf

# find condensation point and change
# VMRs to follow SVP curve 


 svpNa2S = svp_curve(8.550, -13899.0, temp)
 tconNa2S = t_curve(10.045,-0.72,P)
 svpKCl = svp_curve(7.611, -11382.0, temp)
 tconKCl = t_curve(12.479,-0.879,P)

 PPNa = 10**vmr[:,5]*P
 PPK = 10**vmr[:,6]*P

 try:
  vmr[:,3][temp<tconNa2S] = np.log10(svpNa2S[temp<tconNa2S]/P[temp<tconNa2S] )
 except:
  vmr[:,3] = vmr[:,3]

 try:
  vmr[:,4][temp<tconKCl] = np.log10(svpKCl[temp<tconKCl]/P[temp<tconKCl] )
 except:
  vmr[:,4] = vmr[:,4]


 vmr[:,3][vmr[:,3]<-15.0] = -20.0
 vmr[:,4][vmr[:,4]<-15.0] = -20.0
 vmrsum = np.sum(10**vmr[:,:5],axis=1)
 h2 = h2p*(1-vmrsum)
 he = hep*(1-vmrsum)
 vmr[:,5] = np.log10(h2)
 vmr[:,6] = np.log10(he)

# recalculate height profile as VMRs have changed

 (H, temp) = robcat.tempprof(T0,alpha,Teff,tau0,n,vmr,mass,rad)
 if all(H==np.zeros(len(temp)) ):
  return -np.inf

 try:
  Pbase= P[temp<tconNa2S][0]
  Hbase= pchip(P[::-1],H[::-1])
  Hb = np.asscalar(Hbase(Pbase))
  if(Pbase==P[0]):
   return -np.inf
 except:
   return -np.inf

# find folder to calculate in


# run model
 model = nemesisPyMult.nemesispymult(runname, len(spec), len(temp), vmr.shape[1], ith, temp, vmr, mass, rad, H, Hb, opac, FSH, prad, pvar, pimag, preal, Hb2, opac2, FSH2, nav, flat, flon, solzen, emzen, azi, wt)
 
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

FSH2=1.0
Hb2=10.0
opac2=0.0
nav=0
flat=0.0
flon=0.0
solzen=0.0
emzen=0.0
azi=0.0
wt=0.0
pvar=0.1
prad=1.0
pimag=np.array([0.1])
preal=0.1

#print('Reading Na2S Imaginary index data (Norm 1 micron)')
#wvim, imagNa2S = np.loadtxt('Na2S_imag.dat',skiprows=1).T
#wvr, realNa2S = np.loadtxt('Na2S_real.dat',skiprows=1).T

#print('Interpolating imaginary/real data to wavelength range')
#fimag = pchip(wvim,imagNa2S)
#imagwv1 = np.arange(0.8, 1.8, 0.1)
#imagwv2 = np.arange(1.9, 2.9, 0.2)
#imagwv = np.append(imagwv1, imagwv2)
#imagwv = np.append(imagwv,40.0)
#imagwv = np.append(0.2,imagwv)
#pimag = fimag(imagwv)
#freal = pchip(wvr,realNa2S)
#wvnorm = 1.0
#preal = freal(wvnorm)


######################## 
# BEGIN FIGURES
########################

# Param names
titles=np.array(['$\log$NH$_3$', '$\log$CH$_4$', '$\log$H$_2$O', '$\log$Na','$\log$K', 'Mass [M$_{jup}$]','Radius [R$_{jup}$]', '$T_0$', 'Alpha', '$T_{eff}$', 'Tau_0', '$n$', 'Dist', '$\log$FSH','$\log$Opt.','$\log$Part_{Rad.}','$\log$P_{Var.}'])

#Get retrieval info

"""
prefix = 'hr87e/chains/GL570D_'
a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
s=a.get_stats()
values = a.get_equal_weighted_posterior()
"""


values=np.loadtxt('chains/GL570D_phys_live.points')


params=values[:, :n_params]
lnprob=values[:, -1]
samples=params
Nsamp=values.shape[0]
Npars=n_params

NN=10
draws=np.random.randint(len(samples),size=NN)
xrand=samples[draws,:]
pickle.dump(xrand, open('xrandV0.pic','wb') )
model_arr=[]

for i in range(NN):
 mass, rad, T0, alpha, Teff, tau0, n, Pbase, FSH, opac = xrand[i,nvmr:nparams]
 wv, spec, yerr = np.loadtxt("./998/" + runname + ".spx", skiprows=4).T
 model1=calc_spec(xrand[i,:], n_params, 998)
 wv, spec, yerr = np.loadtxt("./997/" + runname + ".spx", skiprows=4).T
 model2=calc_spec(xrand[i,:], n_params, 997)
 model = np.append(model1,model2)
 model_arr = np.concatenate([model_arr, model])
 print model
 


rjup = 69911e3
rad_arr=xrand[:,nvmr+1] * rjup
stfb = 5.670373e-8
model_arr=model_arr.reshape(NN,model.shape[0])
summod = np.sum(model_arr, axis=1)*0.02
print summod

Teff = (summod / (4.0 * np.pi * (rad_arr)**2 * stfb))**0.25
pickle.dump(model_arr, open('model_arr.pic','wb') )
pickle.dump(rad_arr, open('rad_arr.pic','wb') )
pickle.dump(Teff, open('teff_arr.pic','wb') )
print(Teff, np.median(Teff),np.mean(Teff))

