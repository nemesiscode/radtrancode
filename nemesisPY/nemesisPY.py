#!/usr/bin/python
#####################################################################################
#####################################################################################
#      	  			        nemesisPY
#####################################################################################
#####################################################################################

# Version of Nemesis for doing retrievals in solar occultation observations

from nemesis import *
import matplotlib.pyplot as plt
import numpy as np
import nemesisf as ns
import time

runname = raw_input('run name: ')
nmaxproc = 64  #Maximum number of processors to use at the same time

######################################################
#    READING INPUT FILES AND SETTING UP VARIABLES
######################################################

start = time.time()

#Reading .inp
ispace,iscat,ilbl,woff,niter,philimit,nspec,ioff,lin = read_inp_nemesis(runname)
gasgiant = True   #This is assumed. Must be generalized
hcorrx = 0.0
lin = 0            #Does not admit previous retrievals (need to be implemented in this program)



#Reading .spx 
fwhm,xlat,xlon,ngeom,nav,nconv,flat,flon,sol_ang,emiss_ang,azi_ang,wgeom,wave,meas,errmeas = read_spx_nemesis(runname)


#Reading .ref
amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr = read_ref_nemesis(runname)

#Read a priori
lpre = 0
nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac,nx,x0,sx,lx = ns.readapriori(runname,lin,lpre,xlat,npro)

#Reading .set
nmu,mu,wtmu,nf,nphi,isol,dist,lowbc,galb,tsurf,layht,nlayer,laytyp,layint = read_set_nemesis(runname)

#Reading .fla
inormal,iray,ih2o,ich4,io3,inh3,iptf,imie,iuv = read_fla_nemesis(runname)

#Reading .cia 
aname,xdnu,npara = read_cia_nemesis(runname)
flagh2p=0
if npara!=0:
    flagh2p=1


#Loading first state vector
xn = np.zeros([nx])
xa = np.zeros([nx])
sa = np.zeros([nx,nx])
xn[0:nx] = x0[0:nx]
xa[0:nx] = x0[0:nx]
sa[0:nx,0:nx] = sx[0:nx,0:nx]

#Calculating measurement error covariance matrix and measurement vector
ny = sum(nconv)
se = np.zeros([ny,ny])
se1 = np.reshape(np.transpose(errmeas[:,:,0]),[ny])
se1 = se1**2.
y = np.reshape(np.transpose(meas[:,:,0]),[ny])
for i in range(ny):
    se[i,i] = se1[i]   #Assumed to be diagonal


#Read forward modelling errors and add them to se (assumed to be in runname.err file)
#rerr = forwarderr(runname,ngeom,nconv,vconv,woff)    !!!!!!!!!!NEEDS TO BE FINISHED IN nemesis.py 
 
#Checking which elements of the state vector must be fixed
ifix = ns.setifix(x0,sx,nvar,varident,varparam,npro)

#Calculating the calculation wavelengths
maxwave = 1000000
vwavetmp = np.zeros([maxwave,ngeom])
vconv = np.zeros([nconv.max(),ngeom])
nwave = np.zeros([ngeom],dtype='int')
for i in range(ngeom):
    nwave1,vwave1 = wavesetb_nemesis(runname,nconv[i],wave[0:nconv[i],i,0],fwhm)
    vconv[0:nconv[i],i] = wave[0:nconv[i],i,0]
    vwavetmp[0:nwave1,i] = vwave1[0:nwave1]
    nwave[i]=nwave1

#Reading .sur file if required
if gasgiant==False:
    nem,vem,emissivity = read_sur_nemesis(runname)
else:
    nem=2
    vem = np.zeros([nem])
    emissivity = np.zeros([nem])
    vem[0]=-100.0
    vem[1]=1.0e7
    emissivity[0:nem]=1.0


#Opening .itr file to write some information in each iteration so that we can keep track of what's going on
if niter>0:
    fitr = open(runname+'.itr','w')
    fitr.write("\t %i \t %i \t %i\n" % (nx,ny,niter))

######################################################
#    CALCULATE FIRST FORWARD MODEL AND JACOBIAN
######################################################

#Calculating Jacobian
print('nemesisPY :: Calculating Jacobian matrix kk')
ny,yn,kk = jacobian_nemesis(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,ngeom,nav,nconv,vconv,fwhm,ispace,gasgiant,\
              layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,flat,flon,wgeom,lin,nvar,varident,varparam,\
              nx,xn,jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwavetmp,nem,vem,emissivity,nmaxproc,MakePlot=True)


#Calculate gain matrix and average kernels
print('nemesisPY :: Calculating gain matrix')
dd,aa = calc_gain_matrix_nemesis(nx,ny,kk,sa,se)




#Calculate initial value of cost function phi
print('nemesisSO :: Calculating cost function')
chisq,phi = calc_phiret_nemesis(ny,y,yn,se,nx,xn,xa,sa)
ophi = phi
print('chisq/ny = '+str(chisq/float(ny)))



#Assess whether retrieval is likely to be OK
ll = assess_nemesis(nx,ny,kk,sa,se)


######################################################
#        RUN RETRIEVAL FOR EACH ITERATION
######################################################

#Initializing some variables
alambda = 1.0   #Marquardt-Levenberg-type 'braking parameter'
xn1 = np.zeros([nx])
xn1[0:nx] = xn[0:nx]
yn1 = np.zeros([ny])
yn1[0:ny] = yn[0:ny]
for it in range(niter):

    #Writing into .itr file
    fitr.write('%10.5f %10.5f \n' % (chisq,phi))
    for i in range(nx):fitr.write('%10.5f \n' % (xn1[i]))
    for i in range(nx):fitr.write('%10.5f \n' % (xa[i]))
    for i in range(ny):fitr.write('%10.5e \n' % (y[i]))
    for i in range(ny):fitr.write('%10.5e \n' % (se1[i]))
    for i in range(ny):fitr.write('%10.5e \n' % (yn1[i]))
    for i in range(ny):fitr.write('%10.5e \n' % (yn[i]))
    for i in range(nx):
        for j in range(ny):fitr.write('%10.5e \n' % (kk[j,i]))

    #Calculating next state vector
    print('Calculating next iterated state vector')
    x_out = calcnextxn_nemesis(nx,ny,xa,xn,y,yn,dd,aa)
    #  x_out(nx) is the next iterated value of xn using classical N-L
    #  optimal estimation. However, we want to apply a braking parameter
    #  alambda to stop the new trial vector xn1 being too far from the
    #  last 'best-fit' value xn

    for j in range(nx):
        xn1[j] = xn[j] + (x_out[j]-xn[j])/(1.0+alambda)
        #Check to see if log numbers have gone out of range
        #Need to be implemented

    #Test to see if any vmrs have gone negative.
    xflag = 0
    mx = 401
    xnx1 = np.zeros([mx])
    xnx1[0:nx] = xn1[0:nx]
    ncont,xmap,ierr = ns.subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,xlon,nvar,varident,varparam,nx,xnx1,jpre,flagh2p)

    if (ierr==1):
        alambda = almbda * 10.0
        #Need to implement more stuff


    #Calculate test spectrum using trial state vector xn1. 
    #Put output spectrum into temporary spectrum yn1 with
    #temporary kernel matrix kk1. Does it improve the fit? 
#    yn1 = np.zeros([ny])
#    kk1 = np.zeros([ny,nx])
    ny,yn1,kk1 = jacobian_nemesis(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,ngeom,nav,nconv,vconv,fwhm,ispace,gasgiant,\
                                  layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,flat,flon,wgeom,lin,nvar,varident,varparam,\
                                  nx,xn1,jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwavetmp,nem,vem,emissivity,nmaxproc)
 

    #Calculate value of cost function phi
    chisq = 0.0
    chisq,phi = calc_phiret_nemesis(ny,y,yn1,se,nx,xn1,xa,sa)
    print('it.   alambda   ophi   phi')
    print(str(it)+'   '+str(alambda)+'   '+str(ophi)+'   '+str(phi))
    print('chisq/ny = '+str(chisq/float(ny)))
  
    #Does the trial solution fit the data better?
    if (phi <= ophi):
        print('Successful iteration. Updating xn,yn and kk')
        for j in range(nx):
            xn[j] = xn1[j]
        for i in range(ny):
            yn[i] = yn1[i]
            for j in range(nx):
                kk[i,j] = kk1[i,j]


        #Now calculate the gain matrix and averaging kernels
        dd,aa = calc_gain_matrix_nemesis(nx,ny,kk,sa,se)

        #Has solution converged?
        tphi = 100.0*(ophi-phi)/ophi
        if (tphi>=0.0 and tphi<=philimit and alambda<1.0):
            print('phi, phlimit : '+str(tphi)+','+str(philimit))
            print('Phi has converged')
            print('Terminating retrieval')
            break
        else:
           ophi=phi
           alambda = alambda*0.3  #reduce Marquardt brake

    else:
        #Leave xn and kk alone and try again with more braking
        alambda = alambda*10.0  #increase Marquardt brake




#Calculating final covariance matrix
sm,sn,st = calc_serr_nemesis(nx,ny,sa,se,dd,aa)

#Make sure errors stay as a priori for kiter < 0
if niter<0:
    st = sa


######################################################
#            WRITING RESULTS INTO FILES
######################################################


#Closing .itr file
if niter>0:
    fitr.close()


#Writing the .mre 
iform=0
vconv1 = np.transpose(wave[:,:,0])
nspec = 1
ll = write_mre_nemesis(runname,iform,ispace,xlat,xlon,npro,nvar,varident,varparam,nx,ny,y,yn,se,\
                       xa,sa,xn,st,ngeom,nconv,vconv,gasgiant,jpre,jrad,jlogg,iscat,lin,nspec)

#Writing the .cov file
ll = write_cov_nemesis(runname,npro,nvar,varident,varparam,nx,ny,sa,sm,sn,st,se,aa,dd,kk)


#Finishing pogram
end = time.time()
print('Model run OK')
print(' Elapsed time (s) = '+str(end-start))

