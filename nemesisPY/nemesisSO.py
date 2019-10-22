#!/usr/bin/python
#####################################################################################
#####################################################################################
#      	  			        nemesisSO
#####################################################################################
#####################################################################################

# Version of Nemesis for doing retrievals in solar occultation observations

from nemesis import *
import matplotlib.pyplot as plt
import numpy as np
import nemesisf as ns
import time

runname = raw_input('run name: ')
nmaxproc = 32  #Maximum number of processors to use at the same time


#Reading maximum array size in nemesisf in case we need to re-size any array
mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()


######################################################
#    READING INPUT FILES AND SETTING UP VARIABLES
######################################################

start = time.time()

#Reading .inp
ispace,iscat,isolocc,ilbl,inum,ionpeel,woff,niter,philimit,nspec,ioff,lin = read_inp_nemesisl(runname)
gasgiant = False   #This is assumed. Must be generalized
hcorrx = 0.0
lin = 0            #Does not admit previous retrievals (need to be implemented in this program)



#Reading .spx 
fwhm,xlat,xlon,ngeom,nav,wgeom,nconv1,flat,flon,tanhe,wave,meas,errmeas = read_spx_nemesisl(runname)
nconv = nconv1[0]
vconv = np.zeros([nconv])
vconv[:] = wave[:,0,0]   #In solar occultation, it is assumed that all tangent heights have the same convolution wavenumbers




#Reading .ref
amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr = read_ref_nemesis(runname)




#Read a priori
lpre = 0
nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac,nx,x0,sx,lx = ns.readapriori(runname,lin,lpre,xlat,npro)



#Reading .set
nmu,mu,wtmu,nf,nphi,isol,dist,lowbc,galb,tsurf,layht,nlayer,laytyp,layint = read_set_nemesis(runname)
if nlayer != npro-1:
    sys.exit('error in nemesisSO :: Number of layers in .set file must be equal to npro-1 as we are in solar occultations')


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
ny = ngeom*nconv
se = np.zeros([ny,ny])
se1 = np.reshape(np.transpose(errmeas[:,:,0]),[ny])
se1 = se1**2.
y = np.reshape(np.transpose(meas[:,:,0]),[ny])
for i in range(ny):
    se[i,i] = se1[i]   #Assumed to be diagonal



#If options 228 or 229 are set (retrieval of ACS-MIR instrument line shape) then a .fil reference file must be generated
flagils = 0
for i in range(nvar):
    if (varident[i,0] == 228):
        print('nemesisSO :: ILS parameters are to be retrieved using option 228.')
        print('Creating .fil file')
        flagils = 2
        ireq = checkvar_nemesisSO(nx,nvar,npro,varident,varparam)
        iils1 = np.where(ireq == 2)
        iils = iils1[0]
        par1 = xn[iils[0]]
        par2 = xn[iils[1]]
        par3 = xn[iils[2]]
        par4 = xn[iils[3]]
        par5 = xn[iils[4]]
        par6 = xn[iils[5]]
        par7 = xn[iils[6]]
        sys.exit('error :: write_fil_acsmir.f function must be written in python first')
        ll = ns.write_fil_acsmir(runname,nconv,vconv,par1,par2,par3,par4,par5,par6,par7)

    if (varident[i,0] == 229):
        print('nemesisSO :: ILS parameters are to be retrieved using option 229.')
        print('Creating .fil file')
        flagils = 3
        ireq = checkvar_nemesisSO(nx,nvar,npro,varident,varparam)
        iils1 = np.where(ireq == 3)
        iils = iils1[0]
        par1 = xn[iils[0]]
        par2 = xn[iils[1]]
        par3 = xn[iils[2]]
        par4 = xn[iils[3]]
        par5 = xn[iils[4]]
        par6 = xn[iils[5]]
        par7 = xn[iils[6]]
        ll = write_fil_acsmir_v2(runname,nconv,vconv,par1,par2,par3,par4,par5,par6,par7)






#Calculating the calculation wavelengths
if ilbl==2:
    nwave,vwave = wavesetc_nemesis(runname,nconv,vconv,fwhm)
elif ilbl==0:
    nwave,vwave = wavesetb_nemesis(runname,nconv,nconv,fwhm)
else:
    sys.exit('error :: ilbl must be either 0 (correlated-k) or 1 (lbl)')





#Opening .itr file to write some information in each iteration so that we can keep track of what's going on
if niter>0:
    fitr = open(runname+'.itr','w')
    fitr.write("\t %i \t %i \t %i\n" % (nx,ny,niter))




#Checking which elements of the state vector must be fixed
ifix = ns.setifix(x0,sx,nvar,varident,varparam,npro)



#Writing height.lay file (altitude of each layer in the atmosphere)
heightlay = np.zeros([nlayer])
heightlay[0:nlayer] = height[0:npro-1]
ll = write_heightlay_nemesis(nlayer,heightlay)



######################################################
#    CALCULATE FIRST FORWARD MODEL AND JACOBIAN
######################################################


#Calculating Jacobian
print('nemesisSO :: Calculating Jacobian matrix kk')
ny,yn,kk = jacobian_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,
                              ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn,ifix,nmaxproc)


#Calculate gain matrix and average kernels
print('nemesisSO :: Calculating gain matrix')
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
    for i in range(ny):fitr.write('%10.5f \n' % (y[i]))
    for i in range(ny):fitr.write('%10.5f \n' % (se1[i]))
    for i in range(ny):fitr.write('%10.5f \n' % (yn1[i]))
    for i in range(ny):fitr.write('%10.5f \n' % (yn[i]))
    for i in range(nx):
        for j in range(ny):fitr.write('%10.5f \n' % (kk[j,i]))

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
    xnx1 = np.zeros([mx])
    xnx1[0:nx] = xn1[0:nx]
    xflag=0
    ncont,xmap,ierr = ns.subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,xlon,nvar,varident,varparam,nx,xnx1,jpre,flagh2p)

    if (ierr==1):
        alambda = almbda * 10.0
        #Need to implement more stuff


    #Calculate reference .fil file if ILS is to be retrieved
    if flagils == 2:
        iils1 = np.where(ireq == 2)
        iils = iils1[0]
        par1 = xn1[iils[0]]
        par2 = xn1[iils[1]]
        par3 = xn1[iils[2]]
        par4 = xn1[iils[3]]
        par5 = xn1[iils[4]]
        par6 = xn1[iils[5]]
        par7 = xn1[iils[6]]
        ll = ns.write_fil_acsmir(runname,nconv,vconv,par1,par2,par3,par4,par5,par6,par7)

    elif flagils ==3:
        iils1 = np.where(ireq == 3)
        iils = iils1[0]
        par1 = xn1[iils[0]]
        par2 = xn1[iils[1]]
        par3 = xn1[iils[2]]
        par4 = xn1[iils[3]]
        par5 = xn1[iils[4]]
        par6 = xn1[iils[5]]
        par7 = xn1[iils[6]]
        ll = write_fil_acsmir_v2(runname,nconv,vconv,par1,par2,par3,par4,par5,par6,par7)


    #Calculate test spectrum using trial state vector xn1. 
    #Put output spectrum into temporary spectrum yn1 with
    #temporary kernel matrix kk1. Does it improve the fit? 
    ny,yn1,kk1 = jacobian_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,
                              ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn1,ifix,nmaxproc)
 

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
#       CALCULATE FORWARD MODEL FOR EACH GAS
######################################################


specret,specretgas = calcgascn_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn,WRITE_GCN=True)


######################################################
#            WRITING RESULTS INTO FILES
######################################################


#Closing .itr file
if niter>0:
    fitr.close()


#Writing the .mre 
if isolocc==2:
    iform=5
elif isolocc==3:
    iform=5
elif isolocc==4:
    iform=6
else:
    iform=0
vconv1 = wave[:,:,0]
nspec = 1
ll = write_mre_nemesis(runname,iform,ispace,xlat,xlon,npro,nvar,varident,varparam,nx,ny,y,yn,se,\
                       xa,sa,xn,st,ngeom,nconv1,vconv1,gasgiant,jpre,jrad,jlogg,iscat,lin,nspec)

#Writing the .cov file
ll = write_cov_nemesis(runname,npro,nvar,varident,varparam,nx,ny,sa,sm,sn,st,se,aa,dd,kk)


#Finishing pogram
end = time.time()
print('Model run OK')
print(' Elapsed time (s) = '+str(end-start))

