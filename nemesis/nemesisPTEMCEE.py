import numpy as np
import emcee
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import pickle
import triangle
import math
import time
import scipy as sp
import nemesisPTEMCEE
from scipy.stats import invgamma
from scipy.interpolate import interp1d
from math import log,log10
from scipy.stats.kde import gaussian_kde
#                     END OF IMPORTS
#######################################################################
#                    BEGINNING OF DEFINITIONS
#######################################################################

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)

# Define the probability function as likelihood * prior
twopi = 2*math.pi

def lnpriormass(x):

    mass = x[ntemp+nvmr+ngamma]

    if 1 < mass/mjup < 80:
        return 0   #log(1) probabilities are equal
    else:
        return -np.inf  # log(0)  

def lnpriorrad(x):

    rad = x[ntemp+nvmr+ngamma+nmass]

    if 0.6 < rad/rjup < 1.4:
        return 0
    else:
        return -np.inf  # log(0)  


def lnpriorvmr(x):
# vmr = log(vmr)
    vmr = x[xrange(ntemp, ntemp+nvmr)]
    vmrsum = np.sum(10**vmr)
    h2 = h2p*(1-vmrsum)
    he = hep*(1-vmrsum)
    if all(-12<=i<=-1 for i in vmr) and vmrsum <= 0.1:
        return 0   #log(1) probabilities of vmrs are equal (flat?)
    else:
        return -np.inf  # log(0)  

def lnpriorgamma(x):
    gamma=x[ntemp+nvmr]
    invg=invgamma.pdf(gamma,1,scale=5e-5)
    invgprob=invg
    if invgprob > 0:
        return log(invgprob)
    else:
        return -np.inf

def lnpriorb(x):

    b = x[ntemp+nvmr+ngamma+nmass+nrad]
    if 0.01*(min(yerr))**2 <= 10**(b) <= 100*(max(yerr))**2:
        return 0
    else:
        return -np.inf


def lnpriortemp(x):
    temp=x[xrange(ntemp)]
    gamma=x[ntemp+nvmr]
    dt = 0
    for i in range(len(temp)-1):
        if i > 0:
         dt += (temp[i+1]-2*temp[i] + temp[i-1])**2


    lnprobt= -0.5*( dt/gamma + np.log(twopi*gamma) )

    tcubic = cubicinterp(pres0,temp,nempres)

    if np.isfinite(lnprobt) and all(10 <= i <= 7000 for i in tcubic):
        return lnprobt
    else:
        return -np.inf

def lnprior(x):

    lp=lnpriorvmr(x)+lnpriorgamma(x)+lnpriortemp(x)+lnpriorb(x)+lnpriorrad(x)+lnpriormass(x)
#    print(lnpriorvmr(x),lnpriorgamma(x),lnpriortemp(x),lnpriorb(x),lnpriorrad(x),lnpriormass(x))
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp

def cubicinterp(pres,temp,presnew):

# need to have a monotonically increasing x (pres) for function
# so we reverse order of arrays
    pres = pres[::-1]
    temp = temp[::-1]
    f = interp1d(pres, temp, kind='cubic')
    tcubic = f(presnew)
# reverse again to have initial order
    return tcubic

def lnlike(x):

    b = x[ntemp+nvmr+ngamma+nmass+nrad]
    p=np.loadtxt('position.txt')
    ith = -1
    for j in range(len(p)):
     if all(p[j] == x):
      ith = j
      break
    sigma2 = yerr**2 + 10**b
# note we return a single zero
# array when flags prevent VMR etc
# because Fortran will not accept
# a 0th-size array

    if(vflag==0):
     vmr=np.zeros(1)
    else:
     vmr = x[xrange(ntemp, ntemp+nvmr)]
     vmrsum = np.sum(10**vmr)
     h2 = h2p*(1-vmrsum)
     he = hep*(1-vmrsum)
     vmr = np.append(vmr,[np.log10(h2),np.log10(he)])
    if(gflag==0):
     mass=np.zeros(1)
     rad=np.zeros(1)
    else:
     rad = x[ntemp+nvmr+ngamma+nmass]
     mass = x[ntemp+nvmr+ngamma]
    if(tflag==0):
     tcubic=nemtemp
    else:
     temp=x[xrange(ntemp)]    
     tcubic = cubicinterp(pres0,temp,nempres)

#    model = nemesisPTEMCEE.nemesisptemcee(runname, len(y), len(tcubic), len(vmr), ith, tcubic, vmr, vflag, gflag, mass, rad)
    model = 1
    return -0.5*( np.sum( (y-model)**2/sigma2 + np.log(twopi*sigma2) ) )

def lnprob(x):

    lp = lnprior(x)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(x)

def save_object(file, object):
    with open(file, 'wb') as f:
        pickle.dump(object, f)

def load_object(file):
    return pickle.load(open(file, "rb"))

def plot_pdf_kde(ndim, names, samples):
 fig, ax = plt.subplots(figsize=(50, 10))
 rows = 2
 for j in range(ndim):
    plt.subplot(rows,int(math.ceil(float(ndim)/rows)),j)
    pmed = np.percentile(samples[:,j],50)
    pp = np.percentile(samples[:,j],84.14)
    pm = np.percentile(samples[:,j],15.87)
#    plt.hist(sampler.flatchain[:,j], 1000, color="b", histtype="step")
    dist_space = np.linspace( min(samples[:,j]), max(samples[:,j]), 10000 )
    kde = gaussian_kde( samples[:,j] )
    plt.plot(dist_space, kde(dist_space)/max(kde(dist_space)))
    print(pmed, pp, pm)
#    plt.fill_betweenx(kde(dist_space), pm, pp)
#    plt.annotate('%d' % pm,  xy=(pm,0), xycoords='data', fontsize='xx-large')
    plt.axvline(x=pm, linestyle='dashed', color='r', label=r'$%2.2f^{+%2.2f}_{-%2.2f}$' % (pmed, pp-pmed, pmed-pm))
    plt.axvline(x=pp, linestyle='dashed', color='r')
    plt.axvline(x=pmed, linestyle='dashed', color='r')
    plt.title( '%s' % names[j])
#    plt.axis('off')
    plt.legend(loc='best',frameon=False, handlelength=0,numpoints=1, fontsize='xx-large')
    plt.set_xlabel=('Normalised PDF')
    plt.set_ylabel=(r'$\log(VMR)$')

 plt.tight_layout()
 plt.savefig('new.pdf',transparent=True)
 plt.close()

def plot_walkers_iter(ndim, names, chain):

 fig, axis = plt.subplots(figsize=(50, 10))
 rows = 2
 for j in range(ndim):
    plt.subplot(rows,int(math.ceil(float(ndim)/rows)),j)
    plt.plot(chain[:,:,j].T, color="k",alpha=0.8)
    plt.set_xlabel=('n')
    plt.set_ylabel=(names[j])
 plt.tight_layout()
 plt.savefig('lasttime.pdf',transparent=True)
 plt.close()

def plot_spectra(percent, samples, runname, y, yerr, vflag, gflag, ntemp, nvmr, ngamma, nmass, pres0, tempp0, nempres, par_ML):

 x = map(lambda v: [v[1], v[2]-v[1], v[1]-v[0]],
                     zip(*np.percentile(samples, [15.87, 50., 84.14], axis=0)))
 pmparams = np.asarray(x).flatten()

 # median, +1 sigma (high), -1 sigma (low)
 med = pmparams[0::3]
 high = pmparams[1::3]
 low =  pmparams[2::3]
 
 if (tflag==1):
  tempmed = cubicinterp(pres0,med[xrange(ntemp)],nempres)
  temphigh = tempmed + cubicinterp(pres0,high[xrange(ntemp)],nempres)
  templow =  tempmed - cubicinterp(pres0,low[xrange(ntemp)],nempres)  
  aptemp = cubicinterp(pres0,tempp0,nempres)
# maximum likelihood params
  mltemp = cubicinterp(pres0,par_ML[xrange(ntemp)],nempres)
  fig, ax = plt.subplots(figsize=(20, 20))
  plt.semilogy(aptemp, 10**(nempres), label='A Priori', linewidth=4, color='k')
  plt.semilogy(tempmed, 10**(nempres), label='Median', linewidth=4, color='b')
  plt.semilogy(nemtemp, 10**(nempres), label='OE', linewidth=4, color='r')
  plt.semilogy(mltemp, 10**(nempres), label='Max Like', linewidth=4, color='g')
  plt.fill_betweenx(10**(nempres), temphigh, templow, facecolor='b', alpha=0.75)  
  plt.gca().invert_yaxis()
  plt.legend()
  plt.savefig(str(percent)+'%'+"tempprof.pdf")
  plt.close()

 else:
  tempmed, temphigh, templow = nemtemp, nemtemp, nemtemp
 
 if (vflag==1): 
  vmrmed = med[xrange(ntemp, ntemp+nvmr)]
  vmrhigh = vmrmed + high[xrange(ntemp, ntemp+nvmr)]
  vmrlow = vmrmed - low[xrange(ntemp, ntemp+nvmr)]
  vmrsummed, vmrsumhigh, vmrsumlow = np.sum(10**vmrmed), np.sum(10**vmrhigh), np.sum(10**vmrlow) 
  medh2, highh2, lowh2 = h2p*(1-vmrsummed), h2p*(1-vmrsumhigh), h2p*(1-vmrsumlow)
  medhe, highhe, lowhe = hep*(1-vmrsummed), hep*(1-vmrsumhigh), hep*(1-vmrsumlow)

  vmrmed = np.append(vmrmed,[np.log10(medh2),np.log10(medhe)])
  vmrhigh = np.append(vmrhigh,[np.log10(highh2),np.log10(highhe)])
  vmrlow = np.append(vmrlow,[np.log10(lowh2),np.log10(lowhe)])

 else:
  vmrmed, vmrhigh, vmrlow = np.zeros(1), np.zeros(1), np.zeros(1)

 if (gflag==1):
  massmed = med[ntemp+nvmr+ngamma]
  masshigh = massmed + high[ntemp+nvmr+ngamma]
  masslow = massmed - low[ntemp+nvmr+ngamma]
  radmed = med[ntemp+nvmr+ngamma+nmass]
  radhigh = radmed + high[ntemp+nvmr+ngamma+nmass]
  radlow = radmed - low[ntemp+nvmr+ngamma+nmass]
 else: 
   massmed, masshigh, masslow  = 0,0,0
   radmed, radhigh, radlow = 0,0,0  

 sigmamed = np.sqrt(yerr**2 + 10**med[-1]) 
 sigmahigh = np.sqrt(yerr**2 + 10**high[-1])
 sigmalow = np.sqrt(yerr**2 + 10**low[-1])

 #    model = nemesisEMCEE.nemesisemcee(runname, len(y), len(tcubic), len(vmr), ith, tcubic, vmr, vflag, gflag, mass, r$
 print(runname)
 medmodel = nemesisPTEMCEE.nemesisptemcee(runname, len(y), len(tempmed), len(vmrmed), 0, tempmed, vmrmed, vflag, gflag, massmed, radmed)
 highmodel = nemesisPTEMCEE.nemesisptemcee(runname, len(y), len(temphigh), len(vmrhigh), 1, temphigh, vmrhigh, vflag, gflag, masshigh, radhigh)
 lowmodel = nemesisPTEMCEE.nemesisptemcee(runname, len(y), len(templow), len(vmrlow), 2, templow, vmrlow, vflag, gflag, masslow, radlow)

 fig, ax = plt.subplots(figsize=(20, 10))

 plt.plot(wv, y, label='Data', color='g')
 plt.fill_between(wv, (y-sigmamed), (y+sigmamed), facecolor='g', alpha=0.7)

 plt.plot(wv, medmodel, label='Med', color='b')
 plt.fill_between(wv, lowmodel, highmodel, facecolor='b', alpha=0.7)
# plt.fill_between(wv, medmodel-medsigma, medmodel+medsigma, facecolor='b')
# plt.plot(wv, mlmodel, label='Max Like')
 plt.legend()
 plt.tight_layout()
# plt.ylim(0.05*min(y-yerr), 2*max(y+yerr))
 plt.set_xlabel=(r'$Wavelength (\mu m)$')
 plt.set_ylabel=('Luminosity (W)')
 plt.savefig(str(percent)+'%'+"spectra.pdf",transparent=True)
 plt.close()

#                     END OF DEFINITIONS
#######################################################################
#                    BEGINNING OF INPUTS
#######################################################################


# enter runname
runname = 'jupiter_hot'

# read spectrum and errors
input = np.loadtxt("./0/" + runname + ".spx", skiprows=4).T
wv = input[0]
y = input[1]
yerr = input[2]
# read in OE result
nemres = np.loadtxt("./0/" + runname + ".ref", skiprows=13).T
nempres = np.log10(nemres[1])
nemtemp = nemres[2]
# VMRs except H2 He
nemvmr = np.loadtxt("./0/" + runname + ".ref", skiprows=13)[0][3:8]

# a priori / starting values

# take every 5th element of a priori from .ref file as
# starting values, must begin and end at the first and last
# pressures in .ref file!!! 
#pres0=nempres[::6]
#tempp0=nemtemp[::6]
pres0 = np.append(nempres[::6],nempres[-1])
tempp0 = np.append(nemtemp[::6],nemtemp[-1])
vmrp0=np.log10(nemvmr)

gamma0 = 1e5

bp0 = log10(0.01*(min(yerr))**2)
# mp0 is * 1e-24 kg. 1898e24 is mass of jup
mjup = 1898.0
mp0 = 30.0*mjup
# rp0 is in km, rjup
rjup = 69911.0
rp0 = 1.0*rjup

# define h2 and he standard VMRs
h2p = 0.85
hep = 0.15

# define number of params
# no need to change nmass, nrad etc as this will happen at
# flag check stage
ntemp, nvmr, ngamma, nmass, nrad, nb  = len(tempp0), len(vmrp0), 1, 1, 1, 1

# flags for using the VMR and grav(mass+radius).
# if = 0, then the values put in for
# VMR and grav are ignored, so that only
# temp profile can be retrieved. 
# tflag = 0 uses the ref temperature
tflag = 1
vflag = 1
gflag = 1
bflag = 1

if (gflag==0):
 nmass = 0
 nrad = 0 

if (tflag==0):
 ntemp = 0

if(vflag==0):
 nvmr = 0

if(bflag==0):
 nb = 0

ndim = ntemp+nvmr+ngamma+nmass+nrad+nb
# number of walkers
nwalkers = ndim*10
# numer of iterations
niter = 10000
#number of burn-in steps
nburn = 0

#initial positions

p0 = [np.zeros((ndim)) for i in xrange(nwalkers)]

for i in xrange(nwalkers):
 t0 = np.random.random(ntemp)+tempp0
 v0 = np.random.random(nvmr)*(-0.1)+vmrp0
 g0 = np.random.random(1)+gamma0
 m0 = np.random.random(1)+mp0
 r0 = np.random.random(1)+rp0
 b0 = -0.01*bp0*np.random.random(1) + bp0
 p0[i] = np.concatenate((t0,v0,g0,m0,r0,b0))
# p0[i] = v0

continuerun = False
# do not change start
start = 0

if(continuerun == True):
 pos = load_object('pos.bin')
 state = load_object('state.bin')
 prob = load_object('prob.bin')
 start = load_object('start.bin')
 intchain = load_object('chain.bin')
 intflatchain = load_object('flatchain.bin')
 intlikeflat = load_object('likeflat.bin')


#                     END OF INPUTS
#######################################################################
#                    BEGINNING OF MCMC
#######################################################################

# decide number of cores to use here
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=1)
print("Running MCMC...")

# the MCMC will incrementally save probabilities etc
# and make some pretty graphs every 'savenum'
# iterations 
savenum = 100
interval = niter / savenum
if(continuerun == False):
 f = open('maximumprob.dat', 'w')

# do not change sumtime
sumtime = 0

for i in xrange(interval):
# if we are continuing a run, then start will be non-zero
# and it will pick up where it left off (hopefully)
 i = i + start
 print('i, start:')
 print(i, start)

 if (i > interval):
  break

 inittime = time.time()

 if (i == 0):
  pos, prob, state = sampler.run_mcmc(p0, savenum)
 elif (i == nburn/savenum):
# remove burn-in period samples
  print("sampler reset")
  sampler.reset()
  pos, prob, state = sampler.run_mcmc(pos, savenum, lnprob0=prob, rstate0=state)
 else:
  pos, prob, state = sampler.run_mcmc(pos, savenum, lnprob0=prob, rstate0=state)
 f = open('maximumprob.dat', 'a')

 if (continuerun==False):
  samples = sampler.flatchain
  chain = sampler.chain
  lnlike_flat = sampler.flatlnprobability
 else:
  samples = np.append(intflatchain,sampler.flatchain,axis=0)
  chain = np.append(intchain,sampler.chain,axis=1)
  lnlike_flat = np.append(intlikeflat, sampler.flatlnprobability)

# get names of params, then plot PDFs
 tnames = np.array([])
 for t in range(ntemp):
  tna = "T" + str(t+1)
  tnames=np.append(tnames,tna)

 gasnames=np.array(['H2O', 'CO2', 'CO', 'CH4', 'NH3'])
 restnames = np.array(['Gamma', 'Mass', 'Radius', 'b'])
 paramnames = np.append(tnames, gasnames)
 paramnames = np.append(paramnames, restnames)
 plot_pdf_kde(ndim,paramnames,samples)

# work out percentage complete for plot names 
 percent = ((i+1)*100 / interval)

# triangle plot

 figure=triangle.corner(samples,
 quantiles=[0.16,0.5,0.84])

 figure.savefig(str(percent)+ '%'+"triangle.pdf")
 plt.close()

# plot the walker positions over interations

 plot_walkers_iter(ndim, paramnames, chain)

 autocorr_time = sampler.get_autocorr_time()

# get the maximum likelihood value and parameter values 
 f.write('\n\tMean acceptance: %d \n' % (sp.mean(sampler.acceptance_fraction)*100.))
 ML = sp.unravel_index(lnlike_flat.argmax(), lnlike_flat.shape)
 par_ML = samples[ML]
# calculate the parameter medians and +/- 1 sigma uncertainties
# for a normal distribution 68% of scores between +/- 1 sigma
# this corresponds to percentiles of 16 (-1 sig) and 84 (+1 sig)

 f.write('\n\tAutocorrelation time: \n%s' % autocorr_time)
 f.write('\n\tMaximum likelihood parameters = \n%s\nML index = %s \nML value = %s \n' % (par_ML, ML, lnlike_flat[ML]))
# f.write('Median, Upper, Lower = \n%s\n' % vmr)
 f.write('PERCENTAGE COMPLETE: %s%% \n' % percent)

# plot spectra and temperature profile if in retrieval

 plot_spectra(percent, samples, runname, y, yerr, vflag, gflag, ntemp, nvmr, ngamma, nmass, pres0, tempp0, nempres, par_ML)


# save run details to load later
 save_object('state.bin',state)
 save_object('pos.bin', pos)
 save_object('prob.bin', prob)
 save_object('start.bin', i+start+1)
 save_object('chain.bin', chain)
 save_object('flatchain.bin', samples)
 save_object('likeflat.bin', lnlike_flat)
# Calculate Gelman-Rubin statistic  (should be distributed < 1.1 for convergence)
# Take last 50% of chain.
# W is the mean of the variances of each chain.
# chainit is iterations of chain 
 chainlen = len(chain[0,:,0])/2
 wchain = chain[:,chainlen:,:]

 W = np.mean(np.var(wchain, axis=1),axis=0)

# B is the variance of the chain means multipied by n because each chain is based on n draws
# (length of chain = 2n, therefore n is 50% of chain length)

# B = 0.0
# for i in range(nwalkers):
#  B = B + ( np.mean(chain[i,:,:]) - np.mean(chain,axis=0) )**2

 B = (chainlen / (nwalkers-1.0)) * np.var(np.mean(wchain,axis=1),axis=0)
# Variance of stationary distribution is weighted average
# of W and B:

 var = (1.0-(1.0/chainlen))*W + (1.0/chainlen)*B

# because of overdispersion of the starting values,
# this overestimates the true variance, but it is unbiased if
# the starting distribution equals the stationary distribution
# (if starting values were not overdispersed)

# The potential scale reduction factor is
 R = np.sqrt(var/W)

# When R is high (>1.1 or 1.2) then we should run our chains out
# longer to improve convergence to the stationary distribution

 f.write('GELMAN-RUBIN STATISTIC: %s\n' % R)

 totaltime = (time.time() - inittime) / (60*60)
 f.write('TIME TAKEN (HOURS): %s \n' % totaltime)
 sumtime += totaltime
 averagetime = (100.0/percent * sumtime)
 f.write('AVERAGE TIME TO GO (HOURS): %s \n' % averagetime)
 averagetime = averagetime/24.0
 f.write('AVERAGE TIME TO GO (DAYS): %s \n' % averagetime)
 f.close()

print("Done.")

#                     END OF MCMC
#######################################################################

