import numpy as np
import pymultinest
import robcat
import nemesisPyMult
import os
import pdb
import pickle
from scipy.interpolate import pchip
from mpi4py import MPI

######################## 
# CREATE OUTPUT FOLDERS
# GET MPI INFO
########################

if not os.path.exists("chains"): os.mkdir("chains")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

###########################
# DEFINE PRIOR + LIKELIHOOD
###########################


def prior(cube, ndim, nparams):
 
# NH3,CH4, H2O, CO, CO2, Na, K, mass, rad, T0, alpha, Teff, tau0, n, log(Pbase), log(FSH), log(opt.)

# VMRs
 cube[0] = -10.0+9.0*cube[0]
 cube[1] = -10.0+9.0*cube[1]
 cube[2] = -10.0+9.0*cube[2]
 cube[3] = -10.0+9.0*cube[3]
 cube[4] = -10.0+9.0*cube[4]
 cube[5] = -10.0+9.0*cube[5]
 cube[6] = -10.0+9.0*cube[6]

# mass, rad
 cube[7] = 1.0+79.0*cube[7]
 cube[8] = 0.6+1.4*cube[8]

# temperature profile
 cube[9] = 1000.0+7500.0*cube[9]
 cube[10] = 0.5 + 0.5*cube[10]
 cube[11] = 100.0+2000.0*cube[11]
 cube[12] = 7.0*cube[12]
 cube[13] = 1.0+4.0*cube[13]

# gray cloud, log(Pbase), log(FSH), log(opt.)

 cube[14] = np.min(pres) + ( np.max(pres) - np.min(pres) )*cube[14]
 cube[15] = -4 + (1 + 4)*cube[15]
 cube[16] = -6.0 + (4 + 6)*cube[16]

def loglike(cube, ndim, nparams):

# add on H2 and He to VMRs
# ith = rank
# save = np.array(cube[:nparams])
# pickle.dump(save, open( str(ith)+'.pic', 'wb') )
# cube = pickle.load( open('0.pic', 'rb') )
 vmr = np.array(cube[:nvmr])
 vmrsum = np.sum(10**vmr)
 h2 = h2p*(1-vmrsum)
 he = hep*(1-vmrsum)
 vmr = np.append(vmr,[np.log10(h2),np.log10(he)])
 vmr = np.tile(vmr, ( len(pres),1) )

# calculate temperature and height profiles

 mass, rad, T0, alpha, Teff, tau0, n, Pbase, FSH, opac = np.array(cube[nvmr:nparams])
 tau0 = 10**tau0
 FSH = 10**FSH

 (H, temp) = robcat.tempprof(T0,alpha,Teff,tau0,n,vmr,mass,rad)
 if all(H==np.zeros(len(temp)) ):
  return -np.inf

# Find height of cloud

 Hbase= pchip(10**pres[::-1],H[::-1])
 Hb = np.asscalar(Hbase(10**Pbase))

# find folder to calculate in

 ith = rank
# run model
# pdb.set_trace()
 model = nemesisPyMult.nemesispymult(runname, len(spec), len(temp), vmr.shape[1], ith, temp, vmr, mass, rad, H, Hb, opac, FSH, prad, pvar, pimag, preal, Hb2, opac2, FSH2, nav, flat, flon, solzen, emzen, azi, wt)
# calc chi-squared
# pdb.set_trace()
 loglikelihood= -0.5*( np.sum( (spec-model)**2/yerr**2 ) )
 return loglikelihood

######################## 
# BEGIN INPUTS
########################

# enter runname
runname = 'gl570d'

# read spectrum and errors
wv, spec, yerr = np.loadtxt("./0/" + runname + ".spx", skiprows=4).T

# define constants, no. gases, h2 percentage, no. params
nvmr=7

h2p = 0.85
hep = 1.0-h2p

n_params=17

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
# RUN PROGRAM
########################
pymultinest.run(loglike, prior, n_params, outputfiles_basename='./chains/GL570D_',n_live_points=400)

