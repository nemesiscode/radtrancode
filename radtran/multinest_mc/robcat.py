import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import pchip

# NH3,CH4, H2O, CO, CO2, Na, K, H2, He 
k = 1.38e-23
amu = 1.67e-27
R=8.31
G = 6.67408e-11

mmwamu = np.array([17.03052, 16.04246, 18.01528, 28.01, 44.01, 22.9898, 39.0983, 2.01588, 4.002602])

gamma=1.4
P0 = 1000.0
D = 1.5
sigma= 5.670367e-8
mjup = 1898.0e24
rjup = 69911.0e3

def tempprof(T0,alpha,Teff,tau0,n,vmr,mass,rad):

 P = np.logspace(-4,np.log10(P0),80)

 
 beta = alpha*(gamma-1.0)/gamma

# convective

 Tc = T0*(P/P0)**beta

# radiative

 tau = tau0*(P/P0)**n

 Trad = ( Teff**4/2.0 * (1+D*tau) )**0.25

# check to see if realistic, i.e. is the radiative profile stronger than
# convective at the bottom of atmosphere (which should be convective)

 if(Trad[-1]>=Tc[-1]):
  H, Tc, Trad, T = np.zeros(len(P)), np.zeros(len(P)), np.zeros(len(P)),  np.zeros(len(P))
  return (H,T)

# both 

 P = np.log10(P)
 Ptran=np.concatenate( [ P[Trad>Tc][:-5], P[Tc>Trad][5:] ])
 Ttran=np.concatenate( [ Trad[Trad>Tc][:-5], Tc[Tc>Trad][5:] ] )

 Tchip = pchip(Ptran, Ttran)
 T = Tchip(P)

 P = 10**P
# reverse for NEMESIS

 T = T[::-1]
 P = P[::-1]

# calculate height profile
 mmw = np.zeros( len(P) )

 for i in range(len(P)):
  mmw[i] = np.sum(mmwamu * 10**vmr[i]) * amu

 H = np.zeros(len(P))
 SH = np.zeros(len(P))

 Hdiff=np.zeros(len(P)) + 1e6   
 while not all(Hdiff < 0.01):
#constant gravity
 
  radh = rad*rjup + H
  g = G*(mass*mjup)/(radh)**2

  SCALE = (k*T) /(mmw*g)
 
  for i in range(len(P)):
   if (i != 0):
    SH[i]=0.5*(SCALE[i]+SCALE[i-1])
    H[i] = H[i-1] - SH[i] * np.log(P[i]/P[i-1])

  Hzero = H
# variable gravity

  radh=rad*rjup+H
  g = G*(mass*mjup)/(radh)**2
  SCALE = (k*T) /(mmw*g)

  for i in range(len(P)):
   if (i != 0):
    SH[i]=0.5*(SCALE[i]+SCALE[i-1])
    H[i] = H[i-1] - SH[i] * np.log(P[i]/P[i-1])

  Hdiff=abs(H-Hzero)

 H=H/1000.0

 if(any(T<0.0)):
  H, Tc, Trad, T = np.zeros(len(P)), np.zeros(len(P)), np.zeros(len(P)),  np.zeros(len(P))
  return (H,T)

 return (H,T)

#T0,alpha,Teff,tau0,n,vmr,mass,rad=3000.0,0.8,700.0,100.0,1.0,np.array([-3,-3,-3,-4,-4,-4,-4,-1,-0.1]),1.0,1.0
#H,Tc,Trad,T = tempprof(T0,alpha,Teff,tau0,n,vmr,mass,rad)
