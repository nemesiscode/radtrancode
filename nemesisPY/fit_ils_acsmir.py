#!/usr/bin/python
#####################################################################################
#####################################################################################
#                                   fit_ils_acsmir
#####################################################################################
#####################################################################################

#Routine for fitting the ILS in ACS-MIR spectra at one tangent height

from nemesis import *
import matplotlib.pyplot as plt
from shutil import copyfile
import shutil
import numpy as np
import nemesisf as ns
import time

runname = raw_input('run name: ')
alt1 = float(raw_input('Tangent height at which to fit the ILS:'))

nmaxproc = 32  #Maximum number of processors to use at the same time

start = time.time()

######################################################
#        READING THE fit_ils_apr.dat FILE
######################################################

#This file contains the following information:
#  First line: Integer indicating the ILS parametrization (e.g 229)
#  Second line: Number of points for describing parametrization (e.g. 7)
#  Following lines: Each of the parameters and the a priori error

f = open('fit_ils_apr.dat','r')
tmp = np.fromfile(f,sep=' ',count=1,dtype='int')
ilsoption = int(tmp[0])
tmp = np.fromfile(f,sep=' ',count=1,dtype='int')
nils = int(tmp[0])
ilsapr = np.zeros([nils])
ilserr = np.zeros([nils])
for i in range(nils):
    tmp = np.fromfile(f,sep=' ',count=2)
    ilsapr[i] = float(tmp[0])
    ilserr[i] = float(tmp[1])


######################################################
#  CREATING NEW DIRECTORY AND MOVING NECESSARY FILES
######################################################

#Creating directory
path = 'fit_ils'
mkdir_p(path)

#Copying necessary files
copyfile(runname+'.cia', path+'/'+runname+'.cia')
copyfile(runname+'.sur', path+'/'+runname+'.sur')
copyfile(runname+'.ref', path+'/'+runname+'.ref')
copyfile(runname+'.spx', path+'/'+runname+'.spx')
copyfile(runname+'.set', path+'/'+runname+'.set')
copyfile(runname+'.fla', path+'/'+runname+'.fla')
copyfile(runname+'.sha', path+'/'+runname+'.sha')
copyfile(runname+'.xsc', path+'/'+runname+'.xsc')
copyfile(runname+'.lls', path+'/'+runname+'.lls')
copyfile(runname+'.inp', path+'/'+runname+'.inp')
copyfile('aerosol.ref', path+'/'+'aerosol.ref')

#Moving to directory
os.chdir(path)

######################################################
#    READING INPUT FILES AND SETTING UP VARIABLES
######################################################

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

altspx,ialtspx = find_nearest(tanhe,alt1)


#Reading .lls
ngasact,strlta = read_lls_nemesis(runname)


#Reading .ref
amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr = read_ref_nemesis(runname)
altref,ialtref = find_nearest(height,alt1)

#Reading aerosol.ref
npro,naero,height,aerodens = read_aerosol_nemesis()


######################################################
#    WRITING NEW FILES WITH THE REQUIRED SETUP
######################################################


#Writing .ref file
npro = 4
newheight = np.zeros([npro])
newpress = np.zeros([npro])
newtemp = np.zeros([npro])
newvmr = np.zeros([npro,ngas])

ll = write_ref_nemesis(runname,amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,\
                       height[ialtref:ialtref+npro],press[ialtref:ialtref+npro],temp[ialtref:ialtref+npro],vmr[ialtref:ialtref+npro,:])


#Writing aerosol.ref file
ll = write_aerosol_nemesis(npro,naero,height[ialtref:ialtref+npro],aerodens[ialtref:ialtref+npro,:])

#Writing new .spx 
fwhm = -0.1
ngeom = 1
ialtspx1 = int(ialtspx)

ll = write_spx_nemesisSO(runname,fwhm,xlat,xlon,nconv,ngeom,tanhe[ialtspx],wave[:,ialtspx],meas[:,ialtspx],errmeas[:,ialtspx])


#Writing new .inp file (ALL active gases are retrieved + tangent pressure fixed + option 229 for ILS fitting)
fapr = open(runname+'.apr','w')
str1 = 'ILS fit for ACS-MIR'
fapr.write(str1+' \n')

gasactID = np.zeros([ngasact],dtype='int')
isoactID = np.zeros([ngasact],dtype='int')
for i in range(ngasact):
    nwave,vmin,delv,npress,ntemp,gasID1,isoID1,presslevels,templevels = read_ltahead_nemesis(strlta[i])
    gasactID[i] = gasID1
    isoactID[i] = isoID1

nvar = ngasact + 2
fapr.write('\t %i \n' % (nvar))
for i in range(ngasact):
    fapr.write('%i \t %i \t %i \n' % (gasactID[i],isoactID[i],3))
    fapr.write('%f \t %f \n' % (1.0,0.5))

fapr.write('%i \t %i \t %i \n' % (666,0,666))
fapr.write('%7.5f \n' % (altref))
errpress = 1.0e-20
fapr.write('%7.5e \t %7.5e \n' % (press[ialtref],errpress))

fapr.write('%i \t %i \t %i \n' % (ilsoption,0,ilsoption))
for i in range(nils):
    fapr.write('%7.5f \t %7.5f \n' % (ilsapr[i],ilserr[i]))
fapr.close()


######################################################
#    RUNNING THE INVERSE METHOD (COPYING NEMESISSO)
######################################################

#Reading maximum array size in nemesisf in case we need to re-size any array
mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()


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

#Opening .itr file to write some information in each iteration so that we can keep track of what's going on
if niter>0:
    fitr = open(runname+'.itr','w')
    fitr.write("\t %i \t %i \t %i\n" % (nx,ny,niter))


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


#Checking which elements of the state vector must be fixed
ifix = ns.setifix(x0,sx,nvar,varident,varparam,npro)



#Writing height.lay file (altitude of each layer in the atmosphere)
heightlay = np.zeros([nlayer])
heightlay[0:nlayer] = height[0:npro-1]
ll = write_heightlay_nemesis(nlayer,heightlay)


#Calculating Jacobian
print('nemesisSO :: Calculating Jacobian matrix kk')
ny,yn,kk = jacobian_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,\
                              xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn,ifix,nmaxproc)



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
    xflag = 0
    flagh2p = 0
    xnx1 = np.zeros([mx])
    xnx1[0:nx] = xn1[0:nx]
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

    ny,yn1,kk1 = jacobian_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,\
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


#Reading .mre file
lat,lon,ngeom,ny,wave,specret,specmeas,specerrmeas,nx,varident,nxvar,varparam,aprprof,aprerr,retprof,reterr = read_mre_nemesis(runname)


#Writing fit_parameters_ils.dat to store results
os.chdir('../')
fout = open('fit_parameters_ils.dat','w')
str1 = 'Tangent height  ::'
fout.write('%s \t %7.4f \n' % (str1,altref))
str2 = 'Wavenumber range ::'
fout.write('%s \t %7.6f \t %7.6f \n' % (str2,vconv.min(),vconv.max()))
str3 = 'ILS parameters ::'
fout.write(str3+'\n')
for i in range(nils):
    fout.write('\t %7.6f \t  %7.6f \n' % (retprof[i,ngasact+1],reterr[i,ngasact+1]))
fout.close()

#Making plot
axis_font = {'size':'20'}
fig = plt.figure(figsize=(20,8))
wavemin = vconv.min()
wavemax = vconv.max()
ax = plt.axes()
ax.set_xlim(wavemin,wavemax)
ax.tick_params(labelsize=20)
ax.ticklabel_format(useOffset=False)
plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
plt.ylabel('Transmission',**axis_font)
im = ax.plot(wave[:,0],meas[:,0],color='black',linewidth=2.)
ax.fill_between(wave[:,0],specmeas[:,0]-specerrmeas[:,0],specmeas[:,0]+specerrmeas[:,0],alpha=0.3,color='black')
im = ax.plot(vconv[0:nconv],yn[0:nconv],linewidth=2.)
plt.grid()
fig.savefig('ilsfit.png',dpi=100)


#Removing directory
shutil.rmtree(path)


#Finishing pogram
end = time.time()
print('Model run OK')
print(' Elapsed time (s) = '+str(end-start))
