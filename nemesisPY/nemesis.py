#+
# NAME:
#       nemesis.py
#
# DESCRIPTION:
#
#	This library contains functions to read and write files that are formatted as 
#	is the NEMESIS radiative transfer code, as well as functions to replicate some
#       parts of the NEMESIS code
#         
# CATEGORY:
#
#	NEMESIS
# 
# MODIFICATION HISTORY: Juan Alday 21/05/2019

import numpy as np
from struct import * 
import pylab
import sys,os,errno,shutil
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.font_manager as font_manager
import matplotlib as mpl


###############################################################################################

#Defining array sizes that are needed to reconcile between the Fortran and Python routines
#These parameters must be equal to the ones found in arrdef.f and arraylen.f, when the library
#nemesisf was compiled with f2py

def check_arraysize_nemesis():

    """
    FUNCTION NAME : checksize_nemesis()

    DESCRIPTION : Check the maximum array sizes that were used to compile the nemesisf python library
                  with f2py. These parameters must be equal to the ones found in arrdef.f and arraylen.f
                  when the library nemesisf was compiled with f2py.

		  These parameters therefore reflect the sizes of the arrays when any function from the 
                  nemesisf library is called. 

    INPUTS : none 

    OPTIONAL INPUTS: none
            
    OUTPUTS : 

        mx :: Maximum number of parameters in state vector
        my :: Maximum number of parameters in measurement vector
        mconv :: Maximum number of convolution wavenumbers
        mwave :: Maximum number of calculation wavenumbers
        maxpat :: Maximum number of atmospheric paths that can be computed
        maxlay :: Maximum number of layers in the atmosphere
        maxgas :: Maximum number of gases in the atmosphere 
        maxsec :: Maximum number of wavenumbers in aerosol x-section
        mgeom :: Maximum number of geometries
        mvar :: Maximum number of variables to retrieve
        mparam :: Maximum number of extra parameters in varparam 
        maxmu :: Maximum number of zenith angles for scattering calculations

    CALLING SEQUENCE:

        mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    mx = 401
    my = 1024
    mconv = my
    mwave = 20000
    maxout = 500000
    maxout3 = maxout
    maxpat = 60
    maxlay = 2*maxpat
    maxgas = 20
    maxsec = 1000
    mparam = 100
    mvar = 12
    mgeom = 100
    mav = 61
    maxfil = 1000
    mpoint = 16000
    maxmu = 10

    return mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu

###############################################################################################

def find_nearest(array, value):

    """
    FUNCTION NAME : find_nearest()

    DESCRIPTION : Find the closest value in an array

    INPUTS : 

        array :: List of numbers
        value :: Value to search for

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        closest_value :: Closest number to value in array
        index :: Index of closest_value within array

    CALLING SEQUENCE:

        closest_value,index = find_nearest(array,value)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx


###############################################################################################

def read_cia_nemesis(runname):

    """
    FUNCTION NAME : read_cia_nemesis()

    DESCRIPTION : Read the .cia file 

    INPUTS : 

        runname :: Name of the Nemesis run

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
     
        aname :: Name of the CIA file to be used (assumed to exist in the raddata directory)
        xdnu :: Wavenumber step of the CIA table
        npara :: Number of para-H2 fractions listed

    CALLING SEQUENCE:

        aname,xdnu,npara = read_cia_nemesis(runname)

    MODIFICATION HISTORY : Juan Alday (10/06/2019)

    """

    #Opening file
    f = open(runname+'.cia','r')

    s = f.readline().split()   
    aname = s[0]

    s = f.readline().split()
    xdnu = float(s[0])
 
    s = f.readline().split()
    npara = int(s[0])

    return aname,xdnu,npara


###############################################################################################

def read_spx_nemesis(runname, MakePlot=False, SavePlot=False):

    """
    FUNCTION NAME : read_spx_nemesis()

    DESCRIPTION : Reads the .spx file from a Nemesis run

    INPUTS :
 
        runname :: Name of the Nemesis run

    OPTIONAL INPUTS:

        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 
        inst_fwhm :: Instrument full-width at half maximum
        xlat :: Planetocentric latitude at centre of the field of view
        xlon :: Planetocentric longitude at centre of the field of view
        ngeom :: Number of different observation geometries under which the location is observed
        nav(ngeom) ::  For each geometry, nav how many individual spectral calculations need to be
                performed in order to reconstruct the field of fiew
        nconv(ngeom) :: Number of wavenumbers/wavelengths in each spectrum
        flat(ngeom,nav) :: Integration point latitude (when nav > 1)
        flon(ngeom,nav) :: Integration point longitude (when nav > 1)
        sol_ang(ngeom,nav) :: Solar incident angle 
        emiss_ang(ngeom,nav) :: Emission angle
        azi_ang(ngeom,nav) :: Azimuth angle
        wgeom(ngeom,nav) :: Weights for the averaging of the FOV
        wave(nconv,ngeom,nav) :: Wavenumbers/wavelengths
        meas(nconv,ngeom,nav) :: Measured spectrum
        errmeas(nconv,ngeom,nav) :: Measurement noise

    CALLING SEQUENCE:

        inst_fwhm,xlat,xlon,ngeom,nav,nconv,flat,flon,sol_ang,emiss_ang,azi_ang,wgeom,wave,meas,errmeas = read_spx_nemesis(runname)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening file
    f = open(runname+'.spx','r')

    #Reading first line
    tmp = np.fromfile(f,sep=' ',count=4,dtype='float')
    inst_fwhm = float(tmp[0])
    xlat = float(tmp[1])
    xlon = float(tmp[2])
    ngeom = int(tmp[3])

    #Defining variables
    navmax = 1000
    nconvmax = 15000
    nconv = np.zeros([ngeom],dtype='int')
    nav = np.zeros([ngeom],dtype='int')
    flattmp = np.zeros([ngeom,navmax])
    flontmp = np.zeros([ngeom,navmax])
    sol_angtmp = np.zeros([ngeom,navmax])
    emiss_angtmp = np.zeros([ngeom,navmax])
    azi_angtmp = np.zeros([ngeom,navmax])
    wgeomtmp = np.zeros([ngeom,navmax])
    wavetmp = np.zeros([nconvmax,ngeom,navmax])
    meastmp = np.zeros([nconvmax,ngeom,navmax])
    errmeastmp = np.zeros([nconvmax,ngeom,navmax])
    for i in range(ngeom):
        nconv[i] = int(f.readline().strip())
        nav[i] = int(f.readline().strip())
        for j in range(nav[i]):
            tmp = np.fromfile(f,sep=' ',count=6,dtype='float')
            flattmp[i,j] = float(tmp[0])
            flontmp[i,j] = float(tmp[1])
            sol_angtmp[i,j] = float(tmp[2])
            emiss_angtmp[i,j] = float(tmp[3])
            azi_angtmp[i,j] = float(tmp[4])
            wgeomtmp[i,j] = float(tmp[5])
            for iconv in range(nconv[i]):
                tmp = np.fromfile(f,sep=' ',count=3,dtype='float')
                wavetmp[iconv,i,j] = float(tmp[0])
                meastmp[iconv,i,j] = float(tmp[1])
                errmeastmp[iconv,i,j] = float(tmp[2])


    #Making final arrays for the measured spectra
    nconvmax2 = max(nconv)
    navmax2 = max(nav)
    wave = np.zeros([nconvmax2,ngeom,navmax2])
    meas = np.zeros([nconvmax2,ngeom,navmax2])
    errmeas = np.zeros([nconvmax2,ngeom,navmax2])
    flat = np.zeros([ngeom,navmax2])
    flon = np.zeros([ngeom,navmax2])
    sol_ang = np.zeros([ngeom,navmax2])
    emiss_ang = np.zeros([ngeom,navmax2])
    azi_ang = np.zeros([ngeom,navmax2])
    wgeom = np.zeros([ngeom,navmax2])
    for i in range(ngeom):
        wave[0:nconv[i],:,0:nav[i]] = wavetmp[0:nconv[i],:,0:nav[i]]
        meas[0:nconv[i],:,0:nav[i]] = meastmp[0:nconv[i],:,0:nav[i]]
        errmeas[0:nconv[i],:,0:nav[i]] = errmeastmp[0:nconv[i],:,0:nav[i]]  
        flat[:,0:nav[i]] = flattmp[:,0:nav[i]]
        flon[:,0:nav[i]] = flontmp[:,0:nav[i]]
        sol_ang[:,0:nav[i]] = sol_angtmp[:,0:nav[i]]
        emiss_ang[:,0:nav[i]] = emiss_angtmp[:,0:nav[i]]
        azi_ang[:,0:nav[i]] = azi_angtmp[:,0:nav[i]]
        wgeom[:,0:nav[i]] = wgeomtmp[:,0:nav[i]]


    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'size':'13'}
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,figsize=(20,8))
        wavemin = wave.min()
        wavemax = wave.max()
        ax1.set_xlim(wavemin,wavemax)
        ax1.tick_params(labelsize=13)
        ax1.ticklabel_format(useOffset=False)
        ax2.set_xlim(wavemin,wavemax)
        ax2.tick_params(labelsize=13)
        ax2.ticklabel_format(useOffset=False)
        ax2.set_yscale('log')

        ax2.set_xlabel('Wavenumber/Wavelength',**axis_font)
        ax1.set_ylabel('Radiance',**axis_font)  
        ax2.set_ylabel('Radiance',**axis_font)

        for i in range(ngeom):
            im = ax1.plot(wave[0:nconv[i],i,0],meas[0:nconv[i],i,0])
            ax1.fill_between(wave[0:nconv[i],i,0],meas[0:nconv[i],i,0]-errmeas[0:nconv[i],i,0],meas[0:nconv[i],i,0]+errmeas[0:nconv[i],i,0],alpha=0.4)

        for i in range(ngeom):
            im = ax2.plot(wave[0:nconv[i],i,0],meas[0:nconv[i],i,0]) 
            ax2.fill_between(wave[0:nconv[i],i,0],meas[0:nconv[i],i,0]-errmeas[0:nconv[i],i,0],meas[0:nconv[i],i,0]+errmeas[0:nconv[i],i,0],alpha=0.4)
        
        plt.grid()
        plt.show();
        if SavePlot == True:
            fig.savefig(runname+'_spectra.eps',dpi=100)


    return inst_fwhm,xlat,xlon,ngeom,nav,nconv,flat,flon,sol_ang,emiss_ang,azi_ang,wgeom,wave,meas,errmeas



###############################################################################################

def read_spx_nemesisl(runname, MakePlot=False, SavePlot=False):

    """
    FUNCTION NAME : read_spx_nemesisl()

    DESCRIPTION : Reads the .spx file from a NemesisL run

    INPUTS : 
	runname :: Name of the Nemesis run

    OPTIONAL INPUTS:
        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 
        inst_fwhm :: Instrument full-width at half maximum
        xlat :: Planetocentric latitude at centre of the field of view
        xlon :: Planetocentric longitude at centre of the field of view
        ngeom :: Number of different observation geometries under which the location is observed
        nav ::  For each geometry, nav how many individual spectral calculations need to be
                performed in order to reconstruct the field of fiew
        wgeom(ngeom,nav) :: Integration weight
        nconv(ngeom) :: Number of wavenumbers/wavelengths in each spectrum
        flat(ngeom,nav) :: Integration point latitude (when nav > 1)
        flon(ngeom,nav) :: Integration point longitude (when nav > 1)
        tanhe(ngeom) :: Tangent height (km)
        wave(nconv,ngeom,nav) :: Wavenumbers/wavelengths
        meas(nconv,ngeom,nav) :: Measured spectrum
        errmeas(nconv,ngeom,nav) :: Measurement noise

    CALLING SEQUENCE:

	inst_fwhm,xlat,xlon,ngeom,nav,wgeom,nconv,flat,flon,tanhe,wave,meas,errmeas = read_spx_nemesisl(runname)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening file
    f = open(runname+'.spx','r')

    #Reading first line
    tmp = np.fromfile(f,sep=' ',count=4,dtype='float')
    inst_fwhm = float(tmp[0])
    xlat = float(tmp[1])
    xlon = float(tmp[2])
    ngeom = int(tmp[3])

    #Defining variables
    nav = 1 #it needs to be generalized to read more than one NAV per observation geometry 
    nconv = np.zeros([ngeom],dtype='int')
    flat = np.zeros([ngeom,nav])
    flon = np.zeros([ngeom,nav])
    tanhe = np.zeros([ngeom])
    wgeom = np.zeros([ngeom,nav])
    nconvmax = 15000
    wavetmp = np.zeros([nconvmax,ngeom,nav])
    meastmp = np.zeros([nconvmax,ngeom,nav])
    errmeastmp = np.zeros([nconvmax,ngeom,nav])
    for i in range(ngeom):
        nconv[i] = int(f.readline().strip())
        for j in range(nav):
            navsel = int(f.readline().strip())
            tmp = np.fromfile(f,sep=' ',count=6,dtype='float')
            flat[i,j] = float(tmp[0])
            flon[i,j] = float(tmp[1])
            tanhe[i] = float(tmp[2])
            wgeom[i,j] = float(tmp[5])
        for iconv in range(nconv[i]):
            tmp = np.fromfile(f,sep=' ',count=3,dtype='float')
            wavetmp[iconv,i,j] = float(tmp[0])
            meastmp[iconv,i,j] = float(tmp[1])
            errmeastmp[iconv,i,j] = float(tmp[2])


    #Making final arrays for the measured spectra
    nconvmax2 = max(nconv)
    wave = np.zeros([nconvmax2,ngeom,nav])
    meas = np.zeros([nconvmax2,ngeom,nav])
    errmeas = np.zeros([nconvmax2,ngeom,nav])
    for i in range(ngeom):
        wave[0:nconv[i],:,:] = wavetmp[0:nconv[i],:,:]
        meas[0:nconv[i],:,:] = meastmp[0:nconv[i],:,:]
        errmeas[0:nconv[i],:,:] = errmeastmp[0:nconv[i],:,:]

    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'size':'20'}
        fig = plt.figure(figsize=(20,8))
        wavemin = wave.min()
        wavemax = wave.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
#        c = np.linspace(tanhe.min(),tanhe.max(),120)
        nc = 256
        c = np.linspace(0,nc-1,nc)
        tanhelin = np.linspace(tanhe.min(),tanhe.max(),nc)
        norm = mpl.colors.Normalize(vmin=tanhelin.min(), vmax=tanhelin.max())
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
        cmap.set_array([])
        for i in range(ngeom):
            tanhe1,idx = find_nearest(tanhelin,tanhe[i])
            im = ax.plot(wave[0:nconv[i],i,0],meas[0:nconv[i],i,0],c=cmap.to_rgba(idx))
            ax.fill_between(wave[0:nconv[i],i,0],meas[0:nconv[i],i,0]-errmeas[0:nconv[i],i,0],meas[0:nconv[i],i,0]+errmeas[0:nconv[i],i,0],\
                            alpha=0.1,color=cmap.to_rgba(idx))
#            im = ax.plot(wave[0:nconv[i],i,0],meas[0:nconv[i],i,0],color='red')
        cb = plt.colorbar(cmap,ticks=tanhelin)
        tick_locator = mpl.ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        plt.grid()
        plt.show();
        if SavePlot == True:
            fig.savefig(runname+'_spectra.eps',dpi=100)

    return inst_fwhm,xlat,xlon,ngeom,nav,wgeom,nconv,flat,flon,tanhe,wave,meas,errmeas

###############################################################################################

def write_unc_nemesis(runname,npro,tempret,tempreterr,htan,ptan,ptanerr):

    """
    FUNCTION NAME : write_unc_nemesis()

    DESCRIPTION : Write the .unc file for calculating the uncertainties of the pressure levels
                  after the retrieval of a temperature profile and pressure at a given tangent
                  height using the Fortran-written program PTresults

    INPUTS :
 
        runname :: Name of the Nemesis run
        npro :: Number of points in the profile
        tempret(npro) :: Retrieved temperature profile (K)
        tempreterr(npro) :: Uncertainty in the retrieved temperature profile (K)
        htan :: Tangent height at which the pressure has been retrieved (km)
        ptan :: Retrieved pressure at htan (atm)
        ptanerr :: Uncertainty in the retrieved pressure (atm)

    OPTIONAL INPUTS: None
            
    OUTPUTS : 

        Nemesis .unc file 

    CALLING SEQUENCE:

        ll = write_unc_nemesis(runname,npro,tempret,tempreterr,htan,ptan,ptanerr)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    fcdr = open(runname+'.unc','w')

    fcdr.write('%i \n' % (npro))
    for i in range(npro):
        fcdr.write('%7.5f \t %7.5 \n' % (tempret[i],tempreterr[i]))

    fcdr.write('%7.5 \t 10.6e \t 10.6e \n' % (htan,ptan,ptanerr))

    dummy = 1
    return dummy


###############################################################################################

def write_cdr_nemesis(runname,dist,fwhm,ispace,ilbl,nwave,vwave,npath,nconv,vconv,nem,vem,emissivity,tsurf):

    """
    FUNCTION NAME : write_cdr_nemesis()

    DESCRIPTION : Writes a file with .cdr extension, which has the required information to run a 
                  CIRSdrv_wavePY simulation, apart from the standard required files (.prf,.pat,.xsc...) 

    INPUTS :
 
        runname :: Name of the Nemesis run
        dist :: Distance from parent star (AU)
        fwhm :: Full-Width-at-Half-Max
        ispace :: Wavenumbers (0) or wavelengths (1)
        ilbl :: Flag indicating whether to use correlated-k (0) or LBL (2)
        nwave :: Number of calculation wavelengths
        vwave(nwave) :: Calculation wavelengths
        npath :: Number of atmospheric paths to be computed
        nconv :: Number of convolution wavelengths
        vconv(nconv) :: Convolution wavelengths
        nem :: Number of points in surface emissivity spectrum
        vem(nem) :: Wavenumbers for describing the surface emissivity spectrum
        emissivity(nem) :: Surface emissivity spectrum
        tsurf :: Surface temperature (K)

    OPTIONAL INPUTS: None
            
    OUTPUTS : 

        Nemesis .cdr file 

    CALLING SEQUENCE:

        ll = write_cdr_nemesis(runname,dist,fwhm,ispace,ilbl,nwave,vwave,npath,nconv,vconv,nem,vem,emissivity,tsurf)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    fcdr = open(runname+'.cdr','w')
 
    fcdr.write('%7.5f \n' % (dist))
    fcdr.write('%7.5f \n' % (fwhm))
    fcdr.write('%i \n' % (ispace))
    fcdr.write('%i \n' % (ilbl))
#    fcdr.write('%i \n' % (nwave))
#    for i in range(nwave):
#        fcdr.write('%13.10f \n' % (vwave[i]))    

    fcdr.write('%i \n' % (npath))
    fcdr.write('%i \n' % (nconv))
    for i in range(nconv):
        fcdr.write('%13.10f \n' % (vconv[i]))

    fcdr.write('%i \n' % (nem))
    for i in range(nem):
        fcdr.write('%13.10f \t %7.5f \n' % (vem[i],emissivity[i]))
    fcdr.write('%7.5f \n' % (tsurf))
    fcdr.close()

    dummy = 1
    return dummy


###############################################################################################

def write_spx_nemesisSO(runname,inst_fwhm,xlat,xlon,nconv,ngeom,tanhe,wave,meas,errmeas,MakePlot=False,SavePlot=False):
 
    """
    FUNCTION NAME : write_spx_nemesisSO()

    DESCRIPTION : Writes the .spx file from a NemesisL run

    INPUTS :
 
        runname :: Name of the Nemesis run
        inst_fwhm :: Instrument full-width at half maximum
        xlat :: Planetocentric latitude at centre of the field of view
        xlon :: Planetocentric longitude at centre of the field of view
        nconv :: Number of convolution numbers (assumed to be the same for every altitude)
        ngeom :: Number of different observation geometries under which the location is observed
        tanhe(ngeom) :: Tangent height (km)
        wave(nconv,ngeom) :: Wavenumbers/wavelengths
        meas(nconv,ngeom) :: Measured spectrum
        errmeas(nconv,ngeom) :: Measurement noise

    OPTIONAL INPUTS:
        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 

        Nemesis .spx ile 

    CALLING SEQUENCE:

        ll = write_spx_nemesisSO(runname,inst_fwhm,xlat,xlon,nconv,ngeom,tanhe,wave,meas,errmeas)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    fspx = open(runname+'.spx','w')
    fspx.write('%7.5f \t %7.5f \t %7.5f \t %i \n' % (inst_fwhm,xlat,xlon,ngeom))

    for i in range(ngeom):
        fspx.write('\t %i \n' % (nconv))
        fspx.write('\t %i \n' % (1))
        dummy1 = -1.0
        dummy2 = 180.0
        dummy3 = 1.0
        if (ngeom)==1:
            tanhe1 = tanhe
        else:
            tanhe1 = tanhe[i]
        fspx.write('\t %7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \t \n' % (xlat,xlon,tanhe1,dummy1,dummy2,dummy3))
        for k in range(nconv):
                fspx.write('\t %10.5f \t %20.7f \t %20.7f \n' % (wave[k,i],meas[k,i],errmeas[k,i]))

    fspx.close()
    dummy = 1
    return dummy

###############################################################################################


def read_ref_nemesis(runname, MakePlot=False, SavePlot=False):

    """
    FUNCTION NAME : read_ref_nemesis()

    DESCRIPTION : Reads the .ref file from a Nemesis run

    INPUTS : 
        runname :: Name of the Nemesis run

    OPTIONAL INPUTS:
        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 

        amform :: if amform =1 then assumed that all VMR sum up 1. 
        nplanet :: Planet ID (Mercury=1, Venus=2, Earth=3, Mars=4...)
        xlat :: Planetocentric latitude
        npro :: Number of points in the profile
        ngas :: Number of gases whose volume mixing ratios are included in the file
        molwt :: Mean molecular weight of the atmosphere in grams
        gasID(ngas) :: HITRAN ID of the gas that need to be included
        isoID(ngas) :: ID Number of the isotopologue to include (0 for all)
        height(npro) :: height profile in km
        press(npro) :: pressure profile in atm 
        temp(npro) :: temperature profiles in K 
        vmr(npro,ngas) :: volume mixing ratio of the different
 
    CALLING SEQUENCE:

	amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr = read_ref_nemesis(runname)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening file
    f = open(runname+'.ref','r')

    #Reading first and second lines
    tmp = np.fromfile(f,sep=' ',count=1,dtype='int')
    amform = int(tmp[0])
    tmp = np.fromfile(f,sep=' ',count=1,dtype='int')
    
    #Reading third line
    tmp = np.fromfile(f,sep=' ',count=5,dtype='float')
    nplanet = int(tmp[0])
    xlat = float(tmp[1])
    npro = int(tmp[2])
    ngas = int(tmp[3])
    molwt = float(tmp[4])

    #Reading gases
    gasID = np.zeros(ngas,dtype='int')
    isoID = np.zeros(ngas,dtype='int')
    for i in range(ngas):
        tmp = np.fromfile(f,sep=' ',count=2,dtype='int')
        gasID[i] = int(tmp[0])
        isoID[i] = int(tmp[1])

    #Reading profiles
    height = np.zeros(npro)
    press = np.zeros(npro)
    temp = np.zeros(npro)
    vmr = np.zeros([npro,ngas])
    s = f.readline().split()
    for i in range(npro):
        tmp = np.fromfile(f,sep=' ',count=ngas+3,dtype='float')
        height[i] = float(tmp[0])
        press[i] = float(tmp[1])
        temp[i] = float(tmp[2])
        for j in range(ngas):
            vmr[i,j] = float(tmp[3+j])



    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True,figsize=(20,8))
        ax1.plot(press,height,'k-',linewidth=2.)
        ax1.set_xlabel('Pressure (atm)',**axis_font)
        ax1.tick_params(labelsize=20)
        ax1.set_xscale('log')
        ax1.set_ylabel('Altitude (km)',**axis_font)
        ax2.plot(temp,height,'k-',linewidth=2.)
        ax2.set_xlabel('Temperature (K)',**axis_font)
        ax2.tick_params(labelsize=20)
        for j in range(ngas):
            strgas = read_gas_nemesis(gasID[j],isoID[j])
            ax3.plot(vmr[:,j],height,label=strgas,linewidth=2.)

        ax3.set_xscale('log')
        ax3.set_xlabel('Volume mixing ratio',**axis_font)
        ax3.tick_params(labelsize=20)
        plt.subplots_adjust(left=0.06,bottom=0.12,right=0.82,top=0.96,wspace=0.16,hspace=0.20)
        legend = ax3.legend(bbox_to_anchor=(1.7, 1.02),fontsize=20)
        plt.show()

    return amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr


###############################################################################################


def write_ref_nemesis(runname,amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr, MakePlot=False, SavePlot=False):

    """
    FUNCTION NAME : write_ref_nemesis()

    DESCRIPTION : Reads the .ref file from a Nemesis run

    INPUTS : 

        runname :: Name of the Nemesis run
        amform :: if amform =1 then assumed that all VMR sum up 1. 
        nplanet :: Planet ID (Mercury=1, Venus=2, Earth=3, Mars=4...)
        xlat :: Planetocentric latitude
        npro :: Number of points in the profile
        ngas :: Number of gases whose volume mixing ratios are included in the file
        molwt :: Mean molecular weight of the atmosphere in grams
        gasID(ngas) :: HITRAN ID of the gas that need to be included
        isoID(ngas) :: ID Number of the isotopologue to include (0 for all)
        height(npro) :: height profile in km
        press(npro) :: pressure profile in atm 
        temp(npro) :: temperature profiles in K 
        vmr(npro,ngas) :: volume mixing ratio of the different

    OPTIONAL INPUTS:
        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 

        Nemesis .ref file
 
    CALLING SEQUENCE:

        ll = write_ref_nemesis(runname,amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    fref = open(runname+'.ref','w')
    fref.write('\t %i \n' % (amform)) 
    nlat = 1    #Would need to be updated to include more latitudes
    fref.write('\t %i \n' % (nlat))
    
    fref.write('\t %i \t %7.4f \t %i \t %i \t %7.4f \n' % (nplanet,xlat,npro,ngas,molwt))

    gasname = [''] * ngas
    header = [''] * (3+ngas)
    header[0] = 'height(km)'
    header[1] = 'press(atm)'
    header[2] = 'temp(K)  '
    str1 = header[0]+'\t'+header[1]+'\t'+header[2]
    for i in range(ngas):
        fref.write('\t %i \t %i\n' % (gasID[i],isoID[i]))
        strgas = 'GAS'+str(i+1)+'_vmr'
        str1 = str1+'\t'+strgas

    fref.write(str1+'\n')


    for i in range(npro):
        str1 = str('{0:7.6f}'.format(height[i]))+'\t'+str('{0:7.6e}'.format(press[i]))+'\t'+str('{0:7.4f}'.format(temp[i]))
        for j in range(ngas):
           str1 = str1+'\t'+str('{0:7.6e}'.format(vmr[i,j])) 
        fref.write(str1+'\n')
 

    fref.close()
    dummy = 1
    return dummy


###############################################################################################


def write_aerosol_nemesis(npro,naero,height,aerodens, MakePlot=False, SavePlot=False):

    """
    FUNCTION NAME : write_aerosol_nemesis()

    DESCRIPTION : Writes the aerosol.ref file from a Nemesis run

    INPUTS : 

       npro :: Number of points in the profile
       naero :: Number of aerosol types
       height(npro) :: Altitude (km)
       aerodens(npro,naero) :: Aerosol density of each particle type (particles per gram of air)

    OPTIONAL INPUTS:

        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 

        Nemesis aerosol.ref file 
 
    CALLING SEQUENCE:

        ll = write_aerosol_nemesis(npro,naero,height,aerodens)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)
    """

    fref = open('aerosol.ref','w')
    str1 = '# aerosol.ref'
    fref.write(str1+' \n')

    fref.write('\t %i \t %i \n' % (npro,naero))

    for i in range(npro):
        str1 = str('{0:7.6f}'.format(height[i]))
        for j in range(naero):
            str1 = str1+'\t'+str('{0:7.6e}'.format(aerodens[i,j]))
        fref.write(str1+'\n')

    fref.close()
    dummy = 1
    return dummy

###############################################################################################


def read_aerosol_nemesis(MakePlot=False, SavePlot=False):

    """
    FUNCTION NAME : read_aerosol_nemesis()

    DESCRIPTION : Reads the aerosol.ref file from a Nemesis run

    INPUTS : none

    OPTIONAL INPUTS:

        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 

       npro :: Number of points in the profile
       naero :: Number of aerosol types
       height(npro) :: Altitude (km)
       aerodens(npro,naero) :: Aerosol density of each particle type (particles per gram of air)
  
    CALLING SEQUENCE:

        npro,naero,height,aerodens = read_aerosol_nemesis()

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening file
    f = open('aerosol.ref','r')

    #Reading header
    s = f.readline().split()

    #Reading first line
    tmp = np.fromfile(f,sep=' ',count=2,dtype='int') 
    npro = tmp[0]
    naero = tmp[1]

    #Reading data
    height = np.zeros([npro])
    aerodens = np.zeros([npro,naero])
    for i in range(npro):
        tmp = np.fromfile(f,sep=' ',count=naero+1,dtype='float')
        height[i] = tmp[0]
        for j in range(naero):
            aerodens[i,j] = tmp[j+1]


    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        fig = plt.figure(figsize=(10,15))
        ax = plt.axes()
        ax.tick_params(labelsize=20)
        plt.xlabel('Aerosol density (part. per gram of air)',**axis_font)
        plt.ylabel('Altitude (km)',**axis_font)
        for i in range(naero):
                im = ax.plot(aerodens[:,i],height,linewidth=2.)
        plt.grid()
        plt.show()
        if SavePlot == True:
                fig.savefig(runname+'_aerosol.eps',dpi=100)


    return npro,naero,height,aerodens


###############################################################################################

def logflag_nemesis(nvar,npro,varident,varparam):

    """
    FUNCTION NAME : logflag_nemesis()

    DESCRIPTION : It reads the variable ID and returns whether the variable is held in log scale or not

    INPUTS : 

        nvar :: Number of variables
        npro :: Number of points in atmospheric profiles
        varident(nvar,3) :: Variable ID
        varparam(nvar,mparam) :: Extra parameters for describing the retrieved variable

    OPTIONAL INPUTS: none
            
    OUTPUTS : 

        logvar(nx) ::  Flag indicating if element in state vector is in log scale (1) or not (0)

    CALLING SEQUENCE:

        logvar = logflag_nemesis(nvar,npro,varident,varparam)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)
 
    """

    #Reading the number of points associated with each element in state vector
    nxvar = npvar_nemesis(nvar,npro,varident,varparam)
    nx = np.sum(nxvar)
    logvar = np.zeros(nx,dtype='int')

    ix = 0
    for i in range(nvar):
#        if nvar == 1:
#            imod = varident[2]
#            ivar=varident[0]
#        else:
        imod = varident[i,2]
        ivar=varident[i,0]

        for j in range(nxvar[i]):
            iflag = 0
            if ivar!=0:       #For temperature elements are kept linear
                if imod==0:   #Continuous profile
                    iflag=1
                elif imod==1 and j==0:  #Knee profile
                    iflag=1
                elif imod==20 and j==0:  #Knee profile
                    iflag=1
                elif imod==4 and j==0:  #Variable knee profile
                    iflag=1
                elif imod==6 and j==0:  #Venus cloud profile
                    iflag=1
                elif imod==7 and j==0:  #extended profile
                    iflag=1            
                elif imod==17 and j==0:  #extended profile
                    iflag=1
                elif imod==18 and j==0:  #extended profile
                    iflag=1
                elif imod==8 and j==0:  #variable knee profile
                    iflag=1
                elif imod==9 and j==0:  #variable knee profile
                    iflag=1
                elif imod==21 and j==0:  #variable knee profile
                    iflag=1
                elif imod==24 and j==0:  #deep profile
                    iflag=1
                elif imod==27 and j==0:  #step profile (deep value)
                    iflag=1
                elif imod==27 and j==1: #step profile (shallow value)
                    iflag=1
                elif imod==16 and j==0: #lapse rate profile
                    iflag=1
                elif imod==19 and j==0: #lapse rate profile
                    iflag=1
                elif imod==25:  #Shortened continuous model
                    iflag=1
                elif imod==28:  #Modify just one element of a profile
                    iflag=1
         
            if ivar==887:  # X-section spectrum
                iflag=1
            if ivar==888:  # Surface albedo spectrum
                iflag=1
            if ivar==889:  # Surface albedo spectrum multiplier
                iflag=1
            if ivar==666:  # Tangent pressure
                iflag=1
            if ivar==444:  # Particle size and ref. index
                iflag=1
            if ivar==445:  # Particle size and ref. index (coated sphere)
                iflag=1
            if ivar==222:  # Larry's cloud model
                iflag=1
            if ivar==223:  # Larry's revised cloud model
                iflag=1
            if ivar==224:  # Larry's revised cloud model with ext UTC
                iflag=1
            if ivar==225:  # Revised cloud model with ext UTC and trunk.
                iflag=1
            if ivar==226:  # Two cloud model
                iflag=1
            if ivar==227:  # Creme Brulee
                iflag=1


            if imod==1 and j==1:  #log fsh - fixed knee
                iflag=1
            if imod==20 and j==1: #log fsh - fixed knee
                iflag=1
            if imod==4 and j==1: #log fsh - fixed knee
                iflag=1
            if imod==4 and j==2: #variable knee profile
                iflag=1
            if imod==6 and j==1: #Venus cloud profile
                iflag=1
            if imod==7 and j==1: #log fsh - extended
                iflag=1
            if imod==17 and j==1: #log fsh - extended
                iflag=1
            if imod==18 and j==1: #log fsh - extended
                iflag=1
            if imod==8 and j==1: #log fsh - var. knee
                iflag=1
            if imod==8 and j==2: #variable knee profile
                iflag=1
            if imod==9 and j==1: #log fsh - var. knee
                iflag=1
            if imod==19 and j==1: #log fsh - var. knee
                iflag=1
            if imod==9 and j==3: #log cwid - var knee
                iflag=1
            if imod==27 and j==2: #knee pressure for step profile
                iflag=1

            if imod==12 or imod ==13:  #Gaussian/Lorentz cloud
                if j==0:
                    iflag=1
                if j==1:
                    iflag=1
                if j==2:
                    iflag=1

            if imod==14 or imod ==15:  #Gaussian/Lorentz cloud
                if j==0:
                    iflag=1
                if j==2:
                    iflag=1

            if imod==16:  #Lapse rate profile
                if j==1:
                    iflag=1
                if j==2:
                    iflag=1
                if j==3:
                    iflag=1

            if imod==24:  #profile between knee and condensation
                if j==1:
                    iflag=1
                if j==2:
                    iflag=1

            if imod==-1:  #Dust continuous profile in particles/cm3
                iflag=1
            if imod==3: #log scaling factor
                iflag=1
            if imod==10: #log scaling factor
                iflag=1
            if imod==11: #log scaling factor
                iflag=1
            if imod==22: #Brown dwarf T-profile
                iflag=1
            if imod==23: # 2 point vmr gradient
                iflag=1
            if imod==26: #2 point vmr gradient (zero vmr for deep atmos)
                iflag=1

            logvar[ix] = iflag
            ix = ix + 1

    return logvar



###############################################################################################


def read_mre_nemesis(runname,MakePlot=False,SavePlot=False):

    """
    FUNCTION NAME : read_mre_nemesis()

    DESCRIPTION : Reads the .mre file from a Nemesis run

    INPUTS : 
        runname :: Name of the Nemesis run

    OPTIONAL INPUTS:
        MakePlot : If True, a summary plot is made
            
    OUTPUTS : 

 	lat :: Latitude (degrees)
	lon :: Longitude (degrees) 
        ngeom :: Number of geometries in the observation
        ny :: Number of points in the measurement vector for each geometry (assuming they all have the same number of points)
        wave(ny,ngeom) :: Wavelength/wavenumber of each point in the measurement vector
        specret(ny,ngeom) :: Retrieved spectrum for each of the geometries
        specmeas(ny,ngeom) :: Measured spectrum for each of the geometries
        specerrmeas(ny,ngeom) :: Error in the measured spectrum for each of the geometries
        nx :: Number of points in the state vector
        varident(nvar,3) :: Retrieved variable ID, as defined in Nemesis manual
        nxvar :: Number of points in the state vector associated with each retrieved variable
        varparam(nvar,5) :: Extra parameters containing information about how to read the retrieved variables
        aprprof(npro,nvar) :: A priori profile for each variable in the state vector
        aprerr(npro,nvar) :: Error in the a priori profile for each variable in the state vector
        retprof(npro,nvar) :: Retrieved profile for each variable in the state vector
        reterr(npro,nvar) :: Error in the retrieved profile for each variable in the state vector

    CALLING SEQUENCE:

        lat,lon,ngeom,ny,wave,specret,specmeas,specerrmeas,nx,varident,nxvar,varparam,aprprof,aprerr,retprof,reterr = read_mre_nemesis(runname)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening .ref file for getting npro
    amform,nplanet,xlat,npro,ngas,molwt,gasID,isoID,height,press,temp,vmr = read_ref_nemesis(runname) 
    
    #Opening file
    f = open(runname+'.mre','r')

    #Reading first three lines
    tmp = np.fromfile(f,sep=' ',count=1,dtype='int')
    s = f.readline().split()
    nspec = int(tmp[0])
    tmp = np.fromfile(f,sep=' ',count=5,dtype='float')
    s = f.readline().split()
    ispec = int(tmp[0])
    ngeom = int(tmp[1])
    ny2 = int(tmp[2])
    ny = ny2 / ngeom
    nx = int(tmp[3])
    tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
    s = f.readline().split()
    lat = float(tmp[0])
    lon = float(tmp[1])
    
    #Reading spectra
    s = f.readline().split()
    s = f.readline().split()
    wave = np.zeros([ny,ngeom])
    specret = np.zeros([ny,ngeom])
    specmeas = np.zeros([ny,ngeom])
    specerrmeas = np.zeros([ny,ngeom])
    for i in range(ngeom):
        for j in range(ny):
            tmp = np.fromfile(f,sep=' ',count=7,dtype='float')
            wave[j,i] = float(tmp[1])
            specret[j,i] = float(tmp[5])
            specmeas[j,i] = float(tmp[2])
            specerrmeas[j,i] = float(tmp[3])

    #Reading the retrieved state vector
    s = f.readline().split()
    nvar = int(s[1])
    nxvar = np.zeros([nvar],dtype='int')
    aprprof1 = np.zeros([nx,nvar])
    aprerr1 = np.zeros([nx,nvar])
    retprof1 = np.zeros([nx,nvar])
    reterr1 = np.zeros([nx,nvar])
    varident = np.zeros([nvar,3],dtype='int')
    varparam = np.zeros([nvar,5])
    for i in range(nvar):
        s = f.readline().split()
        tmp = np.fromfile(f,sep=' ',count=3,dtype='int')
        varident[i,:] = tmp[:]
        tmp = np.fromfile(f,sep=' ',count=5,dtype='float')
        varparam[i,:] = tmp[:]
        s = f.readline().split()
        np1 = npvar_nemesis(1,npro,varident[i,:],varparam[i,:])
        nxvar[i] = np1
        for j in range(nxvar[i]):
            tmp = np.fromfile(f,sep=' ',count=6,dtype='float')
            aprprof1[j,i] = float(tmp[2])
            aprerr1[j,i] = float(tmp[3])
            retprof1[j,i] = float(tmp[4])
            reterr1[j,i] = float(tmp[5])

    aprprof = np.zeros([nxvar.max(),nvar])
    aprerr = np.zeros([nxvar.max(),nvar])
    retprof = np.zeros([nxvar.max(),nvar])
    reterr = np.zeros([nxvar.max(),nvar])

    for i in range(nvar):
        aprprof[0:nxvar[i],i] = aprprof1[0:nxvar[i],i]
        aprerr[0:nxvar[i],i] = aprerr1[0:nxvar[i],i]
        retprof[0:nxvar[i],i] = retprof1[0:nxvar[i],i]
        reterr[0:nxvar[i],i] = reterr1[0:nxvar[i],i]

 
    return lat,lon,ngeom,ny,wave,specret,specmeas,specerrmeas,nx,varident,nxvar,varparam,aprprof,aprerr,retprof,reterr



###############################################################################################


def read_gas_nemesis(gasID,isoID):

    """
    FUNCTION NAME : read_gas_nemesis()

    DESCRIPTION : Reads the RADTRAN ID for a certain gas and returns some information about it

    INPUTS : 

	gasID :: RADTRAN gas ID
	isoID :: RADTRAN isotopologue ID (0 for all isotopes)

    OPTIONAL INPUTS: none
            
    OUTPUTS : 

	strgas :: String containing the name of the gas

    CALLING SEQUENCE:

        strgas = read_gas_nemesis(gasID,isoID)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """


    if gasID == 1:
        if isoID == 0:
            strgas = '$H_2O$'
        elif isoID == 1:
            strgas = '$H_2^{16}O$'
        elif isoID == 2:
            strgas = '$H_2^{18}O$'
        elif isoID == 3:
            strgas = '$H_2^{17}O$'
        elif isoID == 4:
            strgas = '$HD^{16}O$'
        elif isoID == 5:
            strgas = '$HD^{18}O$'
        elif isoID == 6:
            strgas = '$HD^{17}O$'
        elif isoID == 7:
            strgas = '$D_2^{16}O$'
        else:
            sys.exit('error :: Isotopologue not in the list')


    elif gasID == 2:
        if isoID == 0:
            strgas = '$CO_2$'
        elif isoID == 1:
            strgas = '$^{12}C^{16}O_2$'
        elif isoID == 2:
            strgas = '$^{13}C^{16}O_2$'
        elif isoID == 3:
            strgas = '$^{16}O^{12}C^{18}O$'
        elif isoID == 4:
            strgas = '$^{16}O^{12}C^{17}O$'
        elif isoID == 5:
            strgas = '$^{16}O^{13}C^{18}O$'
        elif isoID == 6:
            strgas = '$^{16}O^{13}C^{17}O$'
        elif isoID == 7:
            strgas = '$^{12}C^{18}O_2$'
        elif isoID == 8:
            strgas = '$^{17}O^{12}C^{18}O$'
        elif isoID == 9:
            strgas = '$^{12}C^{17}O_2$'
        elif isoID == 10:
            strgas = '$^{13}C^{18}O_2$'
        elif isoID == 11:
            strgas = '$^{17}O^{13}C^{18}O$'
        elif isoID == 12:
            strgas = '$^{13}C^{17}O_2$'
        else:
            sys.exit('error :: Isotopologue not in the list')

    elif gasID == 3:
        if isoID == 0:
            strgas = '$O_3$'
        else:
            sys.exit('error :: Isotopologue not in the list')


    elif gasID == 4:
        if isoID == 0:
            strgas = '$N_2O$'
        else:
            sys.exit('error :: Isotopologue not in the list')

    elif gasID == 5:
        if isoID == 0:
            strgas = '$CO$'
        elif isoID == 1:
            strgas = '$^{12}C^{16}O$'
        elif isoID == 2:
            strgas = '$^{13}C^{16}O$'
        elif isoID == 3:
            strgas = '$^{12}C^{18}O$'
        elif isoID == 4:
            strgas = '$^{12}C^{17}O$'
        elif isoID == 5:
            strgas = '$^{13}C^{18}O$'
        elif isoID == 6:
            strgas = '$^{13}C^{17}O$'
        else:
            sys.exit('error :: Isotopologue not in the list')

    else:
        sys.exit('error :: gas not in the list')


    return strgas


###############################################################################################


def read_ktahead_nemesis(filename):

    """
    FUNCTION NAME : read_ktahead_nemesis()

    DESCRIPTION : Read the header information in a correlated-k look-up table written with the standard format of Nemesis

    INPUTS : 

	filename :: Name of the file (supposed to have a .kta extension)

    OPTIONAL INPUTS: none
            
    OUTPUTS : 

        nwave :: Number of wavelength points 
        vmin :: Minimum wavelength
        delv :: Spectral sampling 
        npress :: Number of pressure levels
        ntemp :: Number of temperature levels
        gasID :: RADTRAN gas ID 
        isoID :: RADTRAN isotopologue ID
        pressleves(np) :: Pressure levels (atm)
        templeves(np) :: Temperature levels (K)
 
    CALLING SEQUENCE:

        nwave,vmin,delv,fwhm,npress,ntemp,ng,gasID,isoID,g_ord,del_g,presslevels,templevels = read_ktahead_nemesis(filename)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening file
    strlen = len(filename)
    if filename[strlen-3:strlen] == 'kta':
        f = open(filename,'r')
    else:
        f = open(filename+'.kta','r')

    irec0 = np.fromfile(f,dtype='int32',count=1)
    nwave = np.fromfile(f,dtype='int32',count=1)
    vmin = np.fromfile(f,dtype='float32',count=1)
    delv = np.fromfile(f,dtype='float32',count=1)
    fwhm = np.fromfile(f,dtype='float32',count=1)
    npress = int(np.fromfile(f,dtype='int32',count=1))
    ntemp = int(np.fromfile(f,dtype='int32',count=1))
    ng = int(np.fromfile(f,dtype='int32',count=1))
    gasID = int(np.fromfile(f,dtype='int32',count=1))
    isoID = int(np.fromfile(f,dtype='int32',count=1))

    g_ord = np.fromfile(f,dtype='float32',count=ng)
    del_g = np.fromfile(f,dtype='float32',count=ng)

    presslevels = np.fromfile(f,dtype='float32',count=npress)

    N1 = abs(ntemp)
    if ntemp < 0:
        templevels = np.zeros([npress,n1])
        for i in range(npress):
            for j in range(n1):
                templevels[i,j] =  np.fromfile(f,dtype='float32',count=1)
    else:
        templevels = np.fromfile(f,dtype='float32',count=ntemp)

    return nwave,vmin,delv,fwhm,npress,ntemp,ng,gasID,isoID,g_ord,del_g,presslevels,templevels


###############################################################################################


def read_ltahead_nemesis(filename):

    """
    FUNCTION NAME : read_ltahead_nemesis()

    DESCRIPTION : Read the header information in a line-by-line look-up table written with the standard format of Nemesis

    INPUTS : 

	filename :: Name of the file (supposed to have a .lta extension)

    OPTIONAL INPUTS: none
            
    OUTPUTS : 

        nwave :: Number of wavelength points 
        vmin :: Minimum wavelength
        delv :: Spectral sampling 
        npress :: Number of pressure levels
        ntemp :: Number of temperature levels
        gasID :: RADTRAN gas ID 
        isoID :: RADTRAN isotopologue ID
        pressleves(np) :: Pressure levels (atm)
        templeves(np) :: Temperature levels (K)
 
    CALLING SEQUENCE:

        nwave,vmin,delv,npress,ntemp,gasID,isoID,presslevels,templevels = read_ltahead_nemesis(filename)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening file
    strlen = len(filename)
    if filename[strlen-3:strlen] == 'lta':
        f = open(filename,'r')
    else:
        f = open(filename+'.lta','r')

    irec0 = np.fromfile(f,dtype='int32',count=1)
    nwave = np.fromfile(f,dtype='int32',count=1)
    vmin = np.fromfile(f,dtype='float32',count=1)
    delv = np.fromfile(f,dtype='float32',count=1)
    npress = int(np.fromfile(f,dtype='int32',count=1))
    ntemp = int(np.fromfile(f,dtype='int32',count=1))
    gasID = int(np.fromfile(f,dtype='int32',count=1))
    isoID = int(np.fromfile(f,dtype='int32',count=1))

    presslevels = np.fromfile(f,dtype='float32',count=npress)
    templevels = np.fromfile(f,dtype='float32',count=ntemp)

    return nwave,vmin,delv,npress,ntemp,gasID,isoID,presslevels,templevels


###############################################################################################


def read_gcn_nemesisSO(runname,MakePlot=False,SavePlot=False,ShowPlot=True):

    """
    FUNCTION NAME : read_gcn_nemesisSO()

    DESCRIPTION : Read the .gcn file for a NemesisSO run, which is expected to have the contribution
		  from each gas to the total transmission

    INPUTS : 

	runname :: Name of the Nemesis run

    OPTIONAL INPUTS: 

	MakePlot :: If True, then it creates one plot for each tangent height
        SavePlot :: If True, it saves the figures into .eps files
        ShowPlot :: If False, then figures are not shown in a window, but are kept in the background
            
    OUTPUTS : 

        ny :: Number of convolution wavelengths in each acquisition
        ngeom :: Number of acquisitions in the observation
        ngasact :: Number of active gases that are included in the atmosphere
        gasID(ngasact) :: HITRAN ID of the gas that need to be included
        isoID(ngasact) :: ID Number of the isotopologue to include (0 for all)
        tanhe(ngeom) :: Tangent height of each individual acquisition (km)
        wave(ny,ngeom) :: Wavenumber (cm-1)
        specmeas(ny,ngeom) :: Measured spectra 
        specerrmeas(ny,ngeom) :: Measurement uncertainty
        specret_tot(ny,ngeom) :: Modelled spectra including ALL the active gases
        specret_gas(ny,ngeom,ngasact) :: Spectra with the contribution from each for the active gases
 
    CALLING SEQUENCE:

	ny,ngeom,ngasact,gasID,isoID,tanhe,wave,specmeas,specerrmeas,specret_tot,specret_gas = read_gcn_nemesisSO(runname)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """


    #Opening file
    f = open(runname+'.gcn','r')

    #Reading first line
    tmp = np.fromfile(f,sep=' ',count=3,dtype='int')
    ny = tmp[0]
    ngeom = tmp[1]
    ngasact = tmp[2]

    #Readding active gases
    gasID = np.zeros([ngasact],dtype='int')
    isoID = np.zeros([ngasact],dtype='int')
    for j in range(ngasact):
        tmp = np.fromfile(f,sep=' ',count=2,dtype='int')
        gasID[j] = tmp[0]
        isoID[j] = tmp[1]

    #Reading spectra
    tanhe = np.zeros([ngeom])
    wave = np.zeros([ny,ngeom])
    specmeas = np.zeros([ny,ngeom])
    specerrmeas = np.zeros([ny,ngeom])
    specret_tot = np.zeros([ny,ngeom])
    specret_gas = np.zeros([ny,ngeom,ngasact])
    for i in range(ngeom):
        tanhe[i] = np.fromfile(f,sep=' ',count=1,dtype='float')
        for j in range(ny):
            tmp = np.fromfile(f,sep=' ',count=4+ngasact,dtype='float')
            wave[j,i] = tmp[0]
            specmeas[j,i] = tmp[1]
            specerrmeas[j,i] = tmp[2]
            specret_tot[j,i] = tmp[3]
            specret_gas[j,i,:] = tmp[4:4+ngasact]


    if MakePlot == True:
        axis_font = {'size':'12'}
        legend_font = {'size':'12'}
        for i in range(ngeom):
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True,figsize=(15,8))
            wavemin = wave.min()
            wavemax = wave.max()
            ax2.set_xlim(wavemin,wavemax)
            ax1.tick_params(labelsize=12)
            ax2.tick_params(labelsize=12)
            ax3.tick_params(labelsize=12)
            ax3.set_xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
            ax3.ticklabel_format(useOffset=False)
            ax2.set_ylabel('Transmission',**axis_font)
            ax1.set_ylabel('Transmission',**axis_font)
            ax3.set_ylabel('Residuals',**axis_font)
            ax1.set_title('Spectrum '+str(i+1)+' - Tangent height :: '+str(tanhe[i])+' km',**axis_font)
#            im = ax1.plot(wave[:,i],specmeas[:,i],color='black',linewidth=2.,label='Measured')
#            ax1.fill_between(wave[:,i],specmeas[:,i]-specerrmeas[:,i],specmeas[:,i]+specerrmeas[:,i],alpha=0.3,color='black')
            im = ax1.errorbar(wave[:,i],specmeas[:,i],yerr=specerrmeas[:,i],fmt='o',color='black',ms=2.)
            im = ax1.plot(wave[:,i],specret_tot[:,i],color='green',linewidth=2.,label='Modelled')
            legend = ax1.legend(bbox_to_anchor=(1., 1.02),fontsize=14)
            plt.subplots_adjust(left=0.07,bottom=0.10,right=0.87,top=0.92)
            for j in range(ngasact):
                strgas = read_gas_nemesis(gasID[j],isoID[j])
                im = ax2.plot(wave[:,i],specret_gas[:,i,j],linewidth=2.,label=strgas)
            legend = ax2.legend(bbox_to_anchor=(1., 1.02),fontsize=14)

            ax3.fill_between(wave[:,i],-specerrmeas[:,i],specerrmeas[:,i],color='#C0C0C0')
            im = ax3.plot(wave[:,i],specmeas[:,i]-specret_tot[:,i],color='green',linewidth=2.,label='Measured')
            ax1.grid()
            ax2.grid()
            ax3.grid()
            if ShowPlot==True: 
                plt.show()
            if SavePlot == True:
                fig.savefig(runname+'_gcn_spectra'+str(i+1)+'.eps',dpi=100,transparent=True)

    return ny,ngeom,ngasact,gasID,isoID,tanhe,wave,specmeas,specerrmeas,specret_tot,specret_gas


###############################################################################################

def wavesetb_nemesis(runname,nconv,vconv,fwhm):

    """
    FUNCTION NAME : wavesetb_nemesis()

    DESCRIPTION : Subroutine to calculate which 'calculation' wavelengths are needed to cover the required 'convolution wavelengths'.

    INPUTS : 

        runname :: Name of the Nemesis run
        nconv :: Number of convolution wavelengths
        vconv(nconv) :: Convolution wavelengths
        fwhm :: FWHM of convolved spectrum

    OPTIONAL INPUTS:  none

    OUTPUTS : 

	nwave :: Number of calculation wavenumbers
	vwave(mwave) :: Calculation wavenumbers
 
    CALLING SEQUENCE:

	nwave,vwave = wavesetb_nemesis(runname,nconv,vconv,fwhm)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Reading the .kls file to get the initial and end wavenumbers in the .kta files
    ngasact,strlta = read_kls_nemesis(runname)
    nwavelta = np.zeros([ngasact],dtype='int')
    
    for i in range(ngasact):
        nwave,vmin,delv,fwhmk,npress,ntemp,ng,gasID,isoID,g_ord,del_g,presslevels,templevels = read_ktahead_nemesis(strlta[i])
        nwavelta[i] = nwave

    if len(np.unique(nwavelta)) != 1:
        sys.exit('error :: Number of wavenumbers in all .kta files must be the same')

    vkstart = vmin
    vkstep = delv
    vkend = vkstart + delv*(nwave-1)  

    #Determining the calculation numbers
    savemax = 1000000
    save = np.zeros([savemax])
    ico = 0

    if (vkstep < 0.0 or fwhm == 0.0):
        ico = nconv
        for i in range(nconv):
            save[i] = vconv[i]

    if fwhm < 0.0:
        nconv1,vconv1,nfil,vfil,afil = read_fil_nemesis(runname)
        if nconv != nconv1:
            sys.exit('error :: onvolution wavenumbers must be the same in .spx and .fil files')

        for i in range(nconv1):
            vcentral = vconv1[i]
            for j in range(nconv):
                dv = abs(vcentral-vconv[j])
                if dv < 0.00001:
                    j1 = int((vfil[0,i]-vkstart)/vkstep - 1)
                    j2 = int((vfil[nfil[i]-1,i]-vkstart)/vkstep + 1)
                    v1 = vkstart + (j1-1)*vkstep
                    v2 = vkstart + (j2-1)*vkstep
                    if (v1 < vkstart or v2 > vkend):
                        print('warning from wavesetc')
                        print('Channel wavelengths not covered by lbl-tables')
                        print('v1,v2,vkstart,vkend',v1,v2,vkstart,vkend)
                    for k in range(j2-j1):
                        jj = k + j1
                        vj = vkstart + jj*vkstep
                        save[ico]=vj
                        ico = ico + 1


    elif fwhm > 0.0:

        for i in range(nconv):
            j1 = int( (vconv[i]-0.5*fwhm-vkstart)/vkstep )
            j2 = 2 + int( (vconv[i]+0.5*fwhm-vkstart)/vkstep )
            v1 = vkstart + (j1-1)*vkstep
            v2 = vkstart + (j2-1)*vkstep

            if (v1 < vkstart or v2 > vkend):
                print('warning from wavesetc')
                print('Channel wavelengths not covered by lbl-tables')
                print('v1,v2,vkstart,vkend',v1,v2,vkstart,vkend)

            for k in range(j2-j1):
                jj = k + j1
                vj = v1 + (jj-j1)*vkstep
                save[ico]=vj
                ico = ico + 1

    nco = ico
    #sort calculation wavenumbers into order
    save1 = np.zeros([nco])
    save1 = np.sort(save[0:nco])
  
    #creating calculation wavnumber array
    nwave = nco
    vwave = np.zeros([nwave])
    vwave[:] = save1[:]

    #Now weed out repeated wavenumbers
    vwave[1]=save[1]
    xdiff = 0.9*vkstep  
    ico = 0
    for i in range(nco-1):
        test = abs(save1[i+1]-vwave[ico])
        if test >= xdiff:
            ico = ico + 1
            vwave[ico] = save1[i+1]
            nwave = ico

    print('wavesetb_nemesis :: nwave = '+str(nwave))

    return nwave,vwave



###############################################################################################

def wavesetc_nemesis(runname,nconv,vconv,fwhm):


    """
    FUNCTION NAME : wavesetc_nemesis()

    DESCRIPTION : Subroutine to calculate which 'calculation' wavelengths are needed to cover the required 'convolution wavelengths'.

    INPUTS : 

        runname :: Name of the Nemesis run
        nconv :: NUmber of convolution wavelengths
        vconv(nconv) :: Convolution wavelengths
        fwhm :: Full width at half maximum of instrument line shape

    OPTIONAL INPUTS:  none

    OUTPUTS : 

	nwave :: Number of calculation wavenumbers
	vwave(mwave) :: Calculation wavenumbers
 
    CALLING SEQUENCE:

	nwave,vwave = wavesetc_nemesis(runname,nconv,vconv,fwhm)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Reading the .lls file to get the initial and end wavenumbers in the .lta files
    ngasact,strlta = read_lls_nemesis(runname)
    nwavelta = np.zeros([ngasact],dtype='int')
    
    for i in range(ngasact):
        nwave,vmin,delv,npress,ntemp,gasID,isoID,presslevels,templevels = read_ltahead_nemesis(strlta[i])
        nwavelta[i] = nwave

    if len(np.unique(nwavelta)) != 1:
        sys.exit('error :: Number of wavenumbers in all .lta files must be the same')

    vkstart = vmin
    vkstep = delv
    vkend = vkstart + delv*(nwave-1)    

    #Determining the calculation numbers
    savemax = 1000000
    save = np.zeros([savemax])
    ico = 0
    if fwhm < 0.0:
        nconv1,vconv1,nfil,vfil,afil = read_fil_nemesis(runname)
        if nconv != nconv1:
            sys.exit('error :: convolution wavenumbers must be the same in .spx and .fil files')

        for i in range(nconv1):
            vcentral = vconv1[i]
            for j in range(nconv):
                dv = abs(vcentral-vconv[j])
                if dv < 0.00001:
                    j1 = int((vfil[0,i]-vkstart)/vkstep - 1)
                    j2 = int((vfil[nfil[i]-1,i]-vkstart)/vkstep + 1)
                    v1 = vkstart + (j1-1)*vkstep
                    v2 = vkstart + (j2-1)*vkstep
                    if (v1 < vkstart or v2 > vkend):
                        print('warning from wavesetc')
                        print('Channel wavelengths not covered by lbl-tables')
                        print('v1,v2,vkstart,vkend',v1,v2,vkstart,vkend)
                    for k in range(j2-j1):
                        jj = k + j1
                        vj = vkstart + jj*vkstep
                        save[ico]=vj
                        ico = ico + 1


    elif fwhm > 0.0:
        ishape = read_sha_nemesis(runname)
        if ishape == 0:
            dv = 0.5*fwhm
        elif ishape == 1:
            dv = fwhm
        elif ishape == 2:
            dv = 3.* 0.5 * fwhm / np.sqrt(np.log(2.0))
        else:
            dv = 3.*fwhm

        for i in range(nconv):
            j1 = int( (vconv[i]-dv-vkstart)/vkstep )
            j2 = 2 + int( (vconv[i]+dv-vkstart)/vkstep )
            v1 = vkstart + (j1-1)*vkstep
            v2 = vkstart + (j2-1)*vkstep

            if (v1 < vkstart or v2 > vkend):
                print('warning from wavesetc')
                print('Channel wavelengths not covered by lbl-tables')
                print('v1,v2,vkstart,vkend',v1,v2,vkstart,vkend)

            for k in range(j2-j1):
                jj = k + j1
                vj = v1 + (jj-j1)*vkstep
                save[ico]=vj
                ico = ico + 1

 
    nco = ico
    #sort calculation wavenumbers into order
    save1 = np.zeros([nco])
    save1 = np.sort(save[0:nco])
  
    #creating calculation wavnumber array
    nwave = nco
    vwave = np.zeros([nwave])
    vwave[:] = save1[:]

    #Now weed out repeated wavenumbers
    vwave[1]=save[1]
    xdiff = 0.9*vkstep  
    ico = 0
    for i in range(nco-1):
        test = abs(save1[i+1]-vwave[ico])
        if test >= xdiff:
            ico = ico + 1
            vwave[ico] = save1[i+1]
            nwave = ico

    return nwave,vwave


###############################################################################################


def read_fil_nemesis(runname):

    """
    FUNCTION NAME : read_fil_nemesis()

    DESCRIPTION : Read the .fil file and store the data into variables

    INPUTS : 

        runname :: Name of the Nemesis run

    OPTIONAL INPUTS: none

    OUTPUTS : 

        nconv :: Number of convolution wavelengths 
        wave(nconv) :: Wavenumber array of the spectrum (cm-1)
        nfil(nconv) :: Number of wavelengths used to describe the ILS for each
                        spectral point
        vfil(nfil,nconv) :: Wavenumber array used for describing the ILS (cm-1)
        afil(nfil,nconv) :: Function describing the ILS for each spectral point
  
    CALLING SEQUENCE:

        nconv,wave,nfil,vfil,afil = read_fil_nemesis(runname)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Opening file
    f = open(runname+'.fil','r')

    #Reading first and second lines
    nconv = int(np.fromfile(f,sep=' ',count=1,dtype='int'))
    wave = np.zeros([nconv],dtype='d')
    nfil = np.zeros([nconv],dtype='int')
    nfilmax = 100000
    vfil1 = np.zeros([nfilmax,nconv],dtype='d')
    afil1 = np.zeros([nfilmax,nconv],dtype='d')
    for i in range(nconv):
        wave[i] = np.fromfile(f,sep=' ',count=1,dtype='d')
        nfil[i] = np.fromfile(f,sep=' ',count=1,dtype='int')
        for j in range(nfil[i]):
            tmp = np.fromfile(f,sep=' ',count=2,dtype='d')
            vfil1[j,i] = tmp[0]
            afil1[j,i] = tmp[1]

    nfil1 = nfil.max()
    vfil = np.zeros([nfil1,nconv],dtype='d')
    afil = np.zeros([nfil1,nconv],dtype='d')
    for i in range(nconv):
        vfil[0:nfil[i],i] = vfil1[0:nfil[i],i]
        afil[0:nfil[i],i] = afil1[0:nfil[i],i]

    return nconv,wave,nfil,vfil,afil



###############################################################################################


def read_lls_nemesis(runname):

    """
    FUNCTION NAME : read_lls_nemesis()

    DESCRIPTION : Read the .lls file

    INPUTS : 

        runname :: Name of the Nemesis run

    OPTIONAL INPUTS: none

    OUTPUTS : 

	ngasact :: Number of active gases
	strlta(ngasact) :: String containg the .lta lbl-tables    

    CALLING SEQUENCE:

        ngasact,strlta = read_lls_nemesis(runname)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    ngasact = len(open(runname+'.lls').readlines(  ))

    #Opening file
    f = open(runname+'.lls','r')   
    #strlta = np.chararray(ngasact,itemsize=1000)
    strlta = [''] * ngasact
    for i in range(ngasact):
        s = f.readline().split()
        strlta[i] = s[0]

    return ngasact,strlta

###############################################################################################

def read_kls_nemesis(runname):

    """
    FUNCTION NAME : read_kls_nemesis()

    DESCRIPTION : Read the .kls file

    INPUTS : 

        runname :: Name of the Nemesis run

    OPTIONAL INPUTS: none

    OUTPUTS : 

        ngasact :: Number of active gases
        strkta(ngasact) :: String containg the .kta lbl-tables    

    CALLING SEQUENCE:

        ngasact,strkta = read_kls_nemesis(runname)
 
    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    ngasact = len(open(runname+'.kls').readlines(  ))

    #Opening file
    f = open(runname+'.kls','r')
    #strlta = np.chararray(ngasact,itemsize=1000)
    strkta = [''] * ngasact
    for i in range(ngasact):
        s = f.readline().split()
        strkta[i] = s[0]

    return ngasact,strkta

###############################################################################################


def read_sha_nemesis(runname):

    """
        FUNCTION NAME : read_sha_nemesis()
        
        DESCRIPTION : Read the .sha file
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            lineshape :: Instrument lineshape as defined by Nemesis manual
				(0) Square lineshape
				(1) Triangular
				(2) Gaussian
				(3) Hamming
				(4) Hanning       
 
        CALLING SEQUENCE:
        
            lineshape = read_sha_nemesis(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Opening file
    f = open(runname+'.sha','r')
    s = f.readline().split()
    lineshape = int(s[0])

    return lineshape


###############################################################################################


def read_fla_nemesis(runname):
    
    """
        FUNCTION NAME : read_fla_nemesis()
        
        DESCRIPTION : Read the .fla file
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            runname :: Name of the Nemesis run
            inormal :: ortho/para-H2 ratio is in equilibrium (0) or normal 3:1 (1)
            iray :: (0) Rayleigh scattering optical depth not included
                    (1) Rayleigh optical depths for gas giant atmosphere
                    (2) Rayleigh optical depth suitable for CO2-dominated atmosphere
                    (>2) Rayleigh optical depth suitable for a N2-O2 atmosphere
            ih2o :: Additional H2O continuum off (0) or on (1)
            ich4 :: Additional CH4 continuum off (0) or on (1)
            io3 :: Additional O3 continuum off (0) or on (1)
            inh3 :: Additional NH3 continuum off (0) or on (1)
            iptf :: Normal partition function calculation (0) or high-temperature partition
                    function for CH4 for Hot Jupiters
            imie :: Only relevant for scattering calculations. (0) Phase function is computed
                    from the associated Henyey-Greenstein hgphase*.dat files. (1) Phase function
                    computed from the Mie-Theory calculated PHASEN.DAT
            iuv :: Additional flag for including UV cross sections off (0) or on (1)
        
        CALLING SEQUENCE:
        
            inormal,iray,ih2o,ich4,io3,inh3,iptf,imie,iuv = read_fla_nemesis(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
        """
    
    #Opening file
    f = open(runname+'.fla','r')
    s = f.readline().split()
    inormal = int(s[0])
    s = f.readline().split()
    iray = int(s[0])
    s = f.readline().split()
    ih2o = int(s[0])
    s = f.readline().split()
    ich4 = int(s[0])
    s = f.readline().split()
    io3 = int(s[0])
    s = f.readline().split()
    inh3 = int(s[0])
    s = f.readline().split()
    iptf = int(s[0])
    s = f.readline().split()
    imie = int(s[0])
    s = f.readline().split()
    iuv = int(s[0])
   
    return inormal,iray,ih2o,ich4,io3,inh3,iptf,imie,iuv

###############################################################################################


def read_sur_nemesis(runname):
    
    """
        FUNCTION NAME : read_sur_nemesis()
        
        DESCRIPTION : Read the .sur file (surface emissivity spectrum)
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            nem :: Number of spectral points in surface emissitivity spectrum
            vem(nem) :: Wavenumber array (cm-1)
            emissivity(nem) :: Surface emissivity
        
        CALLING SEQUENCE:
        
            nem,vem,emissivity = read_sur_nemesis(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Opening file
    f = open(runname+'.sur','r')
    nem = int(np.fromfile(f,sep=' ',count=1,dtype='int'))

    vem = np.zeros([nem])
    emissivity = np.zeros([nem])
    for i in range(nem):
        tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
        vem[i] = tmp[0]
        emissivity[i] = tmp[1]

    return nem,vem,emissivity



###############################################################################################


def read_inp_nemesis(runname):
    
    """
        FUNCTION NAME : read_inp_nemesis()
        
        DESCRIPTION : Read the .inp file for a Nemesis run

        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
            iscat :: (0) Thermal emission calculation
                    (1) Multiple scattering required
                    (2) Internal scattered radiation field is calculated first (required for limb-
                        scattering calculations)
                    (3) Single scattering plane-parallel atmosphere calculation
                    (4) Single scattering spherical atmosphere calculation
            ilbl :: (0) Pre-tabulated correlated-k calculation
                    (1) Line by line calculation
                    (2) Pre-tabulated line by line calculation
            
            woff :: Wavenumber/wavelength calibration offset error to be added to the synthetic spectra
            niter :: Number of iterations of the retrieval model required
            philimit :: Percentage convergence limit. If the percentage reduction of the cost function phi
                        is less than philimit then the retrieval is deemed to have converged.
            nspec :: Number of retrievals to perform (for measurements contained in the .spx file)
            ioff :: Index of the first spectrum to fit (in case that nspec > 1).
            lin :: Integer indicating whether the results from previous retrievals are to be used to set any
                    of the atmospheric profiles. (Look Nemesis manual)
        
        CALLING SEQUENCE:
        
            ispace,iscat,ilbl,woff,niter,philimit,nspec,ioff,lin = read_inp_nemesis(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Opening file
    f = open(runname+'.inp','r')
    tmp = f.readline().split()
    ispace = int(tmp[0])
    iscat = int(tmp[1])
    ilbl = int(tmp[2])

    tmp = f.readline().split()
    woff = float(tmp[0])
    fmerrname = str(f.readline().split())
    tmp = f.readline().split()
    niter = int(tmp[0])
    tmp = f.readline().split()
    philimit = float(tmp[0])

    tmp = f.readline().split()
    nspec = int(tmp[0]) 
    ioff = int(tmp[1])
 
    tmp = f.readline().split()
    lin = int(tmp[0])

    return  ispace,iscat,ilbl,woff,niter,philimit,nspec,ioff,lin 



###############################################################################################


def read_inp_nemesisl(runname):
    
    """
        FUNCTION NAME : read_inp_nemesisl()
        
        DESCRIPTION : Read the .inp file for a NemesisL run (slightly different than Nemesis)
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            runname :: Name of the Nemesis run
            ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
            iscat :: (0) Thermal emission calculation
                    (1) Multiple scattering required
                    (2) Internal scattered radiation field is calculated first (required for limb-
                        scattering calculations)
                    (3) Single scattering plane-parallel atmosphere calculation
                    (4) Single scattering spherical atmosphere calculation
            isolocc :: (0) Limb-viewing observation
                    (1) Solar occultation observation (radiance units)
                    (2) Solar occultation observation (transmission)
            ilbl :: (0) Pre-tabulated correlated-k calculation
                    (1) Line by line calculation
                    (2) Pre-tabulated line by line calculation
            inum :: (0) Gradients are calculated analytically
                    (1) Gradients are calculated using numerical differentiation
            ionpeel :: (0) The retrieval is assumed to be a normal NemesisL run,using multiple tangent
                        heights and calculating all atmospheric paths
                        (1) The retrieval is assumed to be set up with just one tangent height, for doing
                            onion-peeling retrievals
            woff :: Wavenumber/wavelength calibration offset error to be added to the synthetic spectra
            niter :: Number of iterations of the retrieval model required
            philimit :: Percentage convergence limit. If the percentage reduction of the cost function phi
                        is less than philimit then the retrieval is deemed to have converged.
            nspec :: Number of retrievals to perform (for measurements contained in the .spx file)
            ioff :: Index of the first spectrum to fit (in case that nspec > 1).
            lin :: Integer indicating whether the results from previous retrievals are to be used to set any
                    of the atmospheric profiles. (Look Nemesis manual)
        
        CALLING SEQUENCE:
        
            ispace,iscat,isolocc,ilbl,inum,ionpeel,woff,niter,philimit,nspec,ioff,lin = read_inp_nemesisl(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Opening file
    f = open(runname+'.inp','r')
    tmp = np.fromfile(f,sep=' ',count=6,dtype='int')
    ispace = int(tmp[0])
    iscat = int(tmp[1])
    isolocc = int(tmp[2])
    ilbl = int(tmp[3])
    inum = int(tmp[4])
    ionpeel = int(tmp[5])

    woff = float(np.fromfile(f,sep=' ',count=1,dtype='float'))
    fmerrname = str(f.readline().split())
    niter = int(np.fromfile(f,sep=' ',count=1,dtype='int'))
    philimit = float(np.fromfile(f,sep=' ',count=1,dtype='float'))

    tmp = np.fromfile(f,sep=' ',count=2,dtype='int')
    nspec = int(tmp[0])
    ioff = int(tmp[1])

    lin = int(np.fromfile(f,sep=' ',count=1,dtype='int'))

    return ispace,iscat,isolocc,ilbl,inum,ionpeel,woff,niter,philimit,nspec,ioff,lin


###############################################################################################


def read_set_nemesis(runname):
    
    """
        FUNCTION NAME : read_set_nemesis()
        
        DESCRIPTION : Read the .set file
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            nmu :: Number of zenith ordinates
            mu(nmu) :: cos(zenith) points
            wtmu(nmu) :: Quadrature weights
            nf :: Required number of Fourier components
            nphi :: Number of azimuth angles
            isol :: Sunlight on/off
            dist :: Solar distance (AU)
            lowbc :: Lower boundary condition (0 Thermal - 1 Lambertian)
            galb :: Ground albedo
            tsurf :: Surface temperature (if planet is not gasgiant)
            layht :: Base height of lowest layer
            nlayer :: Number of vertical levels to split the atmosphere into
            laytp :: Flag to indicate how layering is perfomed (radtran)
            layint :: Flag to indicate how layer amounts are calculated (radtran)
        
        CALLING SEQUENCE:
        
            nmu,mu,wtmu,nf,nphi,isol,dist,lowbc,galb,tsurf,layht,nlayer,laytp,layint = read_set_nemesis(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Opening file
    f = open(runname+'.set','r')
    dummy = f.readline().split()
    nmu1 = f.readline().split()
    nmu = int(nmu1[5])
    mu = np.zeros([nmu],dtype='d')
    wtmu = np.zeros([nmu],dtype='d')
    for i in range(nmu):
        tmp = np.fromfile(f,sep=' ',count=2,dtype='d')
        mu[i] = tmp[0]
        wtmu[i] = tmp[1]
    
    dummy = f.readline().split()
    nf = int(dummy[5])
    dummy = f.readline().split()
    nphi = int(dummy[8])
    dummy = f.readline().split()
    isol = int(dummy[5])
    dummy = f.readline().split()
    dist = float(dummy[5])
    dummy = f.readline().split()
    lowbc = int(dummy[6])
    dummy = f.readline().split()
    galb = float(dummy[3])
    dummy = f.readline().split()
    tsurf = float(dummy[3])

    dummy = f.readline().split()

    dummy = f.readline().split()
    layht = float(dummy[8])
    dummy = f.readline().split()
    nlayer = int(dummy[5])
    dummy = f.readline().split()
    laytp = int(dummy[3])
    dummy = f.readline().split()
    layint = int(dummy[3])

    return nmu,mu,wtmu,nf,nphi,isol,dist,lowbc,galb,tsurf,layht,nlayer,laytp,layint


###############################################################################################


def read_idl_cirsdrv(runname,MakePlot=False,SavePlot=False):
    
    """
        FUNCTION NAME : read_idl_cirsdrv()
        
        DESCRIPTION : Read the .idl file from a CIRSdrv_wave run
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            npath :: Number of paths
            nconv(npath) :: Number of convolution wavelengths
            wave(max(nconv),npath) :: Wavenumber array
            specret(max(nconv),npath) :: Modelled spectra
        
        CALLING SEQUENCE:
        
            npath,nconv,wave,specret = read_idl_cirsdrv(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    f = open(runname+'.idl','r')
    dummy = f.readline().split()
    npath = int(dummy[0])
    
    nconv = np.zeros([npath],dtype='int')
    maxconv = 10000
    wave1 = np.zeros([maxconv,npath])
    specret1 = np.zeros([maxconv,npath])
    
    for i in range(npath):
        dummy = f.readline().split()
        nconv[i] = int(dummy[0])
        #dummy lines
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        #reading spectra
        for j in range(nconv[i]):
            tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
            wave1[j,i] = tmp[0]
            specret1[j,i] = tmp[1]


    wave = np.zeros([nconv.max(),npath])
    specret = np.zeros([nconv.max(),npath])
    for i in range(npath):
        wave[0:nconv[i],i] = wave1[0:nconv[i],i]
        specret[0:nconv[i],i] = specret1[0:nconv[i],i]
    
    
    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        wavemin = wave.min()
        wavemax = wave.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
        colors = plt.cm.jet(np.linspace(0,1,npath))
        for i in range(npath):
            im = ax.plot(wave[0:nconv[i],i],specret[0:nconv[i],i],color=colors[i])
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fmspectra.eps',dpi=100)
    

    return npath,nconv,wave,specret


###############################################################################################


def read_out_cirsdrv(runname,MakePlot=False,SavePlot=False):
    
    """
        FUNCTION NAME : read_out_cirsdrv()
        
        DESCRIPTION : Read the .out file from a CIRSdrv_wave run (convolved and not-convolved spectra)
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            npath :: Number of paths
            nwave(npath) :: Number of calculation wavelengths
            nconv(npath) :: Number of convolution wavelengths
            wavecalc(max(nwave),npath) :: Wavenumber array
            waveconv(max(nconv),npath) :: Wavenumber array
            specret_noconv(max(nwave),npath) :: Modelled spectra (non-convolved)
            specret(max(nconv),npath) :: Modelled spectra (convolved)
        
        CALLING SEQUENCE:
        
            npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret = read_out_cirsdrv(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
        """
    
    f = open(runname+'.out','r')
    dummy = f.readline().split()
    npath = int(dummy[0])
    
    nconv = np.zeros([npath],dtype='int')
    nwave = np.zeros([npath],dtype='int')
    maxwave = 30000
    maxconv = 10000
    wave1 = np.zeros([maxconv,npath])
    specret1 = np.zeros([maxconv,npath])
    wave2 = np.zeros([maxwave,npath])
    specret2 = np.zeros([maxwave,npath])
    
    for i in range(npath):
        dummy = f.readline().split()
        nwave[i] = int(dummy[0])
        #dummy lines
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        #reading spectra
        for j in range(nwave[i]):
            tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
            wave2[j,i] = tmp[0]
            specret2[j,i] = tmp[1]

        #reading convolved spectra
        dummy = f.readline().split()
        nconv[i] = int(dummy[0])
        #dummy lines
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        dummy = f.readline().split()
        #reading spectra
        for j in range(nconv[i]):
            tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
            wave1[j,i] = tmp[0]
            specret1[j,i] = tmp[1]

    wavecalc = np.zeros([nwave.max(),npath])
    specret_noconv = np.zeros([nwave.max(),npath])
    for i in range(npath):
        wavecalc[0:nwave[i],i] = wave2[0:nwave[i],i]
        specret_noconv[0:nwave[i],i] = specret2[0:nwave[i],i]
    

    waveconv = np.zeros([nconv.max(),npath])
    specret = np.zeros([nconv.max(),npath])
    for i in range(npath):
        waveconv[0:nconv[i],i] = wave1[0:nconv[i],i]
        specret[0:nconv[i],i] = specret1[0:nconv[i],i]


    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        wavemin = wavecalc.min()
        wavemax = wavecalc.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
        colors = plt.cm.jet(np.linspace(0,1,npath))
        for i in range(npath):
            im = ax.plot(wavecalc[0:nwave[i],i],specret_noconv[0:nwave[i],i],color=colors[i])
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fmnovconv_spectra.eps',dpi=100)


    return npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret

###############################################################################################


def CIRSdrv_wavePY(runname,MakePlot=False,SavePlot=False):

    """
        FUNCTION NAME : CIRSdrv_wavePY()
        
        DESCRIPTION : Run CIRSdrv_wavePY run (forward model)
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            npath :: Number of paths
            nwave(npath) :: Number of calculation wavelengths
            nconv(npath) :: Number of convolution wavelengths
            wavecalc(max(nwave),npath) :: Wavenumber array
            waveconv(max(nconv),npath) :: Wavenumber array
            specret_noconv(max(nwave),npath) :: Modelled spectra (non-convolved)
            specret(max(nconv),npath) :: Modelled spectra (convolved)
        
        CALLING SEQUENCE:
        
            npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret = CIRSdrv_wavePY(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """
    import subprocess as sub
    from time import sleep

    #Checking that runname.nam exists
    ex = os.path.isfile(runname+'.nam')
    if ex == False:
        f = open(runname+'.nam','w')
        f.write(runname)
        f.close()


    #Running CIRSdrv_waveSO
    strproc = 'CIRSdrv_wavePY < '+runname+'.nam > test_'+runname+'.prc'
    cirsdrv = sub.Popen(strproc, shell = True, stdout = sub.PIPE)
    cirsdrv_out = cirsdrv.communicate()

    sleep(1.0)

    npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret = read_out_cirsdrv(runname)

    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        wavemin = wavecalc.min()
        wavemax = wavecalc.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
        colors = plt.cm.jet(np.linspace(0,1,npath))
        for i in range(npath):
            im = ax.plot(wavecalc[0:nwave[i],i],specret_noconv[0:nwave[i],i],color=colors[i])
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fmnovconv_spectra.eps',dpi=100)

        fig = plt.figure(figsize=(20,8))
        wavemin = waveconv.min()
        wavemax = waveconv.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
        colors = plt.cm.jet(np.linspace(0,1,npath))
        for i in range(npath):
            im = ax.plot(waveconv[0:nconv[i],i],specret[0:nconv[i],i],color=colors[i])
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fmspectra.eps',dpi=100)


    return npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret


###############################################################################################


def CIRSdrv_waveSO(runname,MakePlot=False,SavePlot=False):

    """
        FUNCTION NAME : CIRSdrv_waveSO()
        
        DESCRIPTION : Run CIRSdrv_waveSO run (forward model for solar occultations)
        
        INPUTS :
        
            runname :: Name of the Nemesis run
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            npath :: Number of paths
            nwave(npath) :: Number of calculation wavelengths
            nconv(npath) :: Number of convolution wavelengths
            wavecalc(max(nwave),npath) :: Wavenumber array
            waveconv(max(nconv),npath) :: Wavenumber array
            specret_noconv(max(nwave),npath) :: Modelled spectra (non-convolved)
            specret(max(nconv),npath) :: Modelled spectra (convolved)
        
        CALLING SEQUENCE:
        
            npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret = CIRSdrv_waveSO(runname)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    import subprocess as sub
    from time import sleep
 
    #Checking that runname.nam exists
    ex = os.path.isfile(runname+'.nam')
    if ex == False:
        f = open(runname+'.nam','w')
        f.write(runname)
        f.close()


    #Running CIRSdrv_waveSO
    strproc = 'CIRSdrv_waveSO < '+runname+'.nam > test_'+runname+'.prc' 
    cirsdrv = sub.Popen(strproc, shell = True, stdout = sub.PIPE)
    cirsdrv_out = cirsdrv.communicate()

    sleep(1.0)

    npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret = read_out_cirsdrv(runname)

    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        wavemin = wavecalc.min()
        wavemax = wavecalc.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
        colors = plt.cm.jet(np.linspace(0,1,npath))
        for i in range(npath):
            im = ax.plot(wavecalc[0:nwave[i],i],specret_noconv[0:nwave[i],i],color=colors[i])
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fmnovconv_spectra.eps',dpi=100)

        fig = plt.figure(figsize=(20,8))
        wavemin = waveconv.min()
        wavemax = waveconv.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
        colors = plt.cm.jet(np.linspace(0,1,npath))
        for i in range(npath):
            im = ax.plot(waveconv[0:nconv[i],i],specret[0:nconv[i],i],color=colors[i])
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fmspectra.eps',dpi=100)


    return npath,nwave,nconv,wavecalc,waveconv,specret_noconv,specret


###############################################################################################


def gsetrad_nemesis(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,nconv,vconv,fwhm,ispace,gasgiant,\
                     layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,xlat,xlon,lin,nvar,varident,varparam,\
                     nx,xn,jalb,jxsc,jtan,jpre,tsurf):


    """
        FUNCTION NAME : gsetrad_nemesis()
        
        DESCRIPTION : Calls the Fortran-written routine gsetrad.f, in order to create the .pat, .prf, .xsc and .sca files for a 
                      CIRSdrv_wavePY simulation. This function only reconciles the array sizes between NemesisPY and the array sizes 
                      used to compile the Nemesis Fortran library with f2py

        INPUTS :
        
            runname :: Name of the Nemesis run
            iscat :: 0=thermal emission; 1=plane parallel scattering; 2= limb/near-limb scattering; 3= single-scattering calculations
            nmu :: Number of zenith ordinates
            mu(nmu) :: Cos(zenith angles)
            wtmu(nmu) :: Quadrature weights
            isol :: Sunlight on(1)/off(0)
            dist :: Solar distance (AU)
            lowbc :: LOwer boundary condition
            galb :: ground albedo
            nf :: Required number of Fourier components
            nconv :: Number of convolution wavelengths
            vconv(nconv) :: Convolution wavelength array
            fwhm :: FWHM of convolved spectrum
            ispace :: 0=cm-1, 1= wavelengths
            gasgiant :: Indicates if planet is a gas giant
            layht :: Altitude of base layer (if non-limb, this is set to -emiss_ang for limb calculations)
            nlayer :: Number of atmospheric layers
            laytyp :: How layers are separated
            layint :: How layer amounts are integrated
            sol_ang :: Solar zenith angle
            emiss_ang :: Thermal emission angle
            azi_ang :: Azimuth angle
            xlat :: latitude of observation
            lin :: Integer to indicate role of previous retrieval (if any)
            nvar :: Number of retrieved variables
            varident(nvar,3) :: identity of constituent to retrieved and parameterisation
            varparam(nvar,3) :: Additional parameters constraining profile
            nx :: Number of elements in state vector
            xn(nx) :: State vector
            jalb :: Position of surface albedo spectrum in xn
            jxsc :: Position of x-section spectrum in xn
            jtan :: Position of tangent height correction in xn
            jpre :: Position of tangent pressure
            tsurf :: Surface temperature

        OPTIONAL INPUTS: none

        OUTPUTS :
        
            xmap(maxv,maxgas+2+maxcon,maxpro) :: Mapping array relate functional derivatives calculated 
                                                 by CIRSRADG to the state vector elements
                                               
            tsurf :: Updated surface temperature (if using previously retrieved value)

 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
    
    """
 
    import nemesisf as ns
 
    #Reading maximum array sizes in the Fortran routines
    mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()

    #Re-sizing arrays for reconciling Fortran and Python codes
    mu1 = np.zeros([maxmu]) 
    wtmu1 = np.zeros([maxmu])
    mu1[0:nmu] = mu[0:nmu]
    wtmu1[0:nmu] = wtmu[0:nmu]

    vconv1 = np.zeros([mconv])
    vconv1[0:nconv] = vconv[0:nconv]

    varident1 = np.zeros([mvar,3],dtype='int')
    varparam1 = np.zeros([mvar,mparam])
    nparam = len(varparam[0,:])
    varident1[0:nvar,:] = varident[0:nvar,:]
    varparam1[0:nvar,0:nparam] = varparam[0:nvar,0:nparam]   

    xn1 = np.zeros([mx])
    xn1[0:nx] = xn[0:nx] 

    xmap = ns.gsetrad(runname,iscat,nmu,mu1,wtmu1,isol,dist,lowbc,galb,nf,nconv,vconv1,fwhm,ispace,gasgiant,\
                      layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,xlat,xlon,lin,nvar,varident1,varparam1,\
                      nx,xn1,jalb,jxsc,jtan,jpre,tsurf)

    return xmap,tsurf


###############################################################################################

def nemesisSOfm_parallel(ix,runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,xlat,xlon,lin,\
                   ilbl,nvar,varident,varparam,jpre,nx,nxn,xnx):

    """
        FUNCTION NAME : nemesisSOfm_parallel()
        
        DESCRIPTION : Reads some input variables like the state vector, and prepares the required files for a CIRSdrv_wavePY run,
                      with some assumption valid for the solar occultation observations.
        
        INPUTS :
       
            ix :: Index indicating which state vector in xnx must be used 
            runname :: Name of the Nemesis run
            nconv :: Number of convolution wavelengths
            vconv(nconv) :: Wavenumber array (cm-1)
            nwave :: Number of calculation wavelengths
            vwave(nwave) :: Calculation wavenumber array (cm-1)
            npro :: Number of altitude levels in atmosphere
            height(npro) :: Altitude (km)
            ngeom :: Number of tangent heights
            tanhe(ngeom) :: Tangent height (km)
            fwhm :: Full width at half maximum (cm-1) 
            ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
            xlat :: Latitude of observation
            xlon :: Longitude of observation
            lin :: Unit number of previous retrieval (if any)
            ilbl :: Flag indicating whether to use correlated-k (0) or lbl (2)
            nvar :: Number of variables to retrieve
            varident(nvar,3) :: Variable ID
            varparam(nvar,5) :: Other parameters defining the retrieved variables 
            jpre :: Position of tangent pressure (if retrieved)
            nx :: Number of elements in state vector
            nxn :: Number of state vectors in xnx
            xn(nx,nxn) :: Matrix containing all state vectors to be computed

        OPTIONAL INPUTS: none
        
        OUTPUTS :
       
            specret(nconv,ngeom) :: Modelled spectra (convolved)
        
        CALLING SEQUENCE:
         
            specret = nemesisSOfm_parallel(ix,runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,nxn,xnx)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    import nemesisf as ns
    from shutil import copyfile
    from time import sleep

    #Parameters assumed for this version of Nemesis
    iscat = 0 #No scattering
    gasgiant = False #Assumed that it is Mars
    isolocc = 2 #Transmission (not thermal emission)
    tsurf = 100.0 #arbitrary
    dist = 1.5 #arbitrary in transmission
    ionpeel = 1
    layht = 0.0
    laytyp = 5
    layint = 1
    hcorrx = 0.0
    nlayer = npro-1

    #Making array for specifying element in the state vector with special requirements
    ireq = checkvar_nemesisSO(nx,nvar,npro,varident,varparam)

    #Calculating the number oh paths that need to be computed based on 
    mintanhe = tanhe.min()
    maxtanhe = tanhe.max()
    jlevlo1 = np.where(height == mintanhe)
    jlevhi1 = np.where(height == maxtanhe)

    if (len(jlevlo1) > 1):
        sys.exit('error in nemesisSOfm_parallel :: jlevlo has more than one value')
    if (len(jlevhi1) > 1):
        sys.exit('error in nemesisSOfm_parallel :: jlevhi has more than one value')

    jlevlo = int(jlevlo1[0]) + 1  #because of the fortran indexing
    jlevhi = int(jlevhi1[0]) + 1

    #Making new arrays so that they have the same dimensions as in the Fortran routines
    mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()
    vconvnem = np.zeros([mconv])
    vconvnem[0:nconv] = vconv[:]

    #Creating directory to calculate the run and moving necessary files
    path = str(ix)
    mkdir_p(path)

    #Copying necessary files
    copyfile(runname+'.cia', path+'/'+runname+'.cia')
    copyfile(runname+'.ref', path+'/'+runname+'.ref')
    copyfile(runname+'.fla', path+'/'+runname+'.fla')
    copyfile(runname+'.xsc', path+'/'+runname+'.xsc')
    copyfile('aerosol.ref', path+'/'+'aerosol.ref')

    if laytyp==5:   #Base altitude of atmospheric layers must be read from height.lay
        copyfile('height.lay', path+'/'+'height.lay')

    if ilbl==2:   #LBL
        copyfile(runname+'.lls', path+'/'+runname+'.lls')
        if fwhm > 0.0:
            copyfile(runname+'.sha', path+'/'+runname+'.sha')
    if ilbl==0:  #Correlated-k
        copyfile(runname+'.kls', path+'/'+runname+'.kls')
        if fwhm < 0.0:
            copyfile(runname+'.fil', path+'/'+runname+'.fil')


    #Checking if ILS is going to be retrieved
    iflag = 0
    for i in range(nx):
        if ireq[i] == 2:  #Option 228
            iflag = 2
        if ireq[i] == 3:
            iflag = 3     #Option 229

    #Creating new .fil file if necessary
    if fwhm < 0.0:  #Creating .fil file in case it is necessary
        if iflag==0:
            copyfile(runname+'.fil', path+'/'+runname+'.fil')
            os.chdir(path)
        if iflag==2:  #Option 228
            os.chdir(path)
            iils1 = np.where(ireq == 2)
            iils = iils1[0]
            par1 = xn1[iils[0]]
            par2 = xn1[iils[1]]
            par3 = xn1[iils[2]]
            par4 = xn1[iils[3]]
            par5 = xn1[iils[4]]
            par6 = xn1[iils[5]]
            par7 = xn1[iils[6]]
            ll = ns.write_fil_acsmir(runname,nconv,vconvnem,par1,par2,par3,par4,par5,par6,par7)
        if iflag==3:  #Option 229
            os.chdir(path)
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
    else:
        os.chdir(path)

    xn1 = np.zeros([mx])
    xn1[0:nx] = xnx[0:nx,ix]

    #Creating files for communicating with CIRSdrv_wavePY
    xmap = ns.gsetradl(runname,nconv,vconvnem,fwhm,ispace,iscat,gasgiant,layht,nlayer,laytyp,layint,xlat,xlon,lin,hcorrx,nvar,varident,varparam,nx,xn1,jpre,tsurf,isolocc,ionpeel,jlevlo,jlevhi)

    #In solar occultation we do not mind about the surface emissivity
    nem = 2
    vem = np.zeros([nem])
    emissivity = np.zeros([nem])
    vem[0]=-100.0
    vem[1]=1.0e7
    emissivity[0:nem]=1.0

    #Writing .cdr file for communicating with CIRSdrv_wavePY
    npath = ngeom
    ll = write_cdr_nemesis(runname,dist,fwhm,ispace,ilbl,nwave,vwave,npath,nconv,vconv,nem,vem,emissivity,tsurf)

    #Running CIRSdrv_wavePY
    npath,nwave1,nconv1,wavecalc,waveconv,specret_noconv,specret = CIRSdrv_wavePY(runname)
    os.chdir('../')

    #Removing folder
    shutil.rmtree(path)

    return specret


###############################################################################################

def nemesisSOfm(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,xlat,xlon,lin,\
                   ilbl,nvar,varident,varparam,jpre,nx,xn):

    """
        FUNCTION NAME : nemesisSO()
        
        DESCRIPTION : Reads some input variables like the state vector, and prepares the required files for a CIRSdrv_wave2 run.
                      It also reads an index, which indicates if some element of the state vector must be perturbed. Finalle,
                      it calls CIRSdrv_waveSO and waits for the output
        
        INPUTS :
        
            runname :: Name of the Nemesis run
            nconv :: Number of convolution wavelengths
            vconv(nconv) :: Wavenumber array (cm-1)
            nwave :: Number of calculation wavelengths
            vwave(nwave) :: Calculation wavenumber array (cm-1)
            npro :: Number of altitude levels in atmosphere
            height(npro) :: Altitude (km)
            ngeom :: Number of tangent heights
            tanhe(ngeom) :: Tangent height (km)
            fwhm :: Full width at half maximum (cm-1) 
            ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
            xlat :: Latitude of observation
            xlon :: Longitude of observation
            lin :: Unit number of previous retrieval (if any)
            ilbl :: Flag indicating whether to use correlated-k (0) or lbl (2)
            nvar :: Number of variables to retrieve
            varident(nvar,3) :: Variable ID
            varparam(nvar,5) :: Other parameters defining the retrieved variables 
            jpre :: Position of tangent pressure (if retrieved)
            nx :: Number of elements in state vector
            xn(nx) :: State vector 

        OPTIONAL INPUTS: none
        
        OUTPUTS :
       
            specret_noconv(nconv,npath) :: Modelled spectra (no convolved) 
            specret(max(nconv),npath) :: Modelled spectra (convolved)
        
        CALLING SEQUENCE:
        
            specret_noconv,specret = nemesisSOfm(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    import nemesisf as ns
    from shutil import copyfile
    from time import sleep

    #Parameters assumed for this version of Nemesis
    iscat = 0 #No scattering
    gasgiant = False #Assumed that it is Mars
    isolocc = 2 #Transmission (not thermal emission)
    tsurf = 100.0 #arbitrary
    dist = 1.5 #arbitrary in transmission
    ionpeel = 1
    layht = 0.0
    laytyp = 5
    layint = 1
    hcorrx = 0.0
    nlayer = npro-1

    #Making array for specifying element in the state vector with special requirements
    ireq = checkvar_nemesisSO(nx,nvar,npro,varident,varparam)

    #Calculating the number oh paths that need to be computed based on 
    mintanhe = tanhe.min()
    maxtanhe = tanhe.max()
    jlevlo1 = np.where(height == mintanhe)
    jlevhi1 = np.where(height == maxtanhe)

    if (len(jlevlo1) > 1):
        sys.exit('error in nemesisSOfm :: jlevlo has more than one value')
    if (len(jlevhi1) > 1):
        sys.exit('error in nemesisSOfm :: jlevhi has more than one value')

    jlevlo = int(jlevlo1[0]) + 1  #because of the fortran indexing
    jlevhi = int(jlevhi1[0]) + 1

    #Making new arrays so that they have the same dimensions as in the Fortran routines
    mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()
    vconvnem = np.zeros([mconv])
    vconvnem[0:nconv] = vconv[:]

    #Creating directory to calculate the run and moving necessary files
    ix=0
    path = str(ix)
    mkdir_p(path)

    #Copying necessary files
    copyfile(runname+'.cia', path+'/'+runname+'.cia')
    copyfile(runname+'.ref', path+'/'+runname+'.ref')
    copyfile(runname+'.fla', path+'/'+runname+'.fla')
    copyfile(runname+'.xsc', path+'/'+runname+'.xsc')
    copyfile('aerosol.ref', path+'/'+'aerosol.ref')

    if laytyp==5:   #Base altitude of atmospheric layers must be read from height.lay
        copyfile('height.lay', path+'/'+'height.lay')

    if ilbl==2:   #LBL
        copyfile(runname+'.lls', path+'/'+runname+'.lls')
        if fwhm > 0.0:
            copyfile(runname+'.sha', path+'/'+runname+'.sha')
        if fwhm < 0.0:
            copyfile(runname+'.fil', path+'/'+runname+'.fil')
    if ilbl==0:  #Correlated-k
        copyfile(runname+'.kls', path+'/'+runname+'.kls')
        if fwhm < 0.0:
            copyfile(runname+'.fil', path+'/'+runname+'.fil')
      

    #Checking if ILS is going to be retrieved
    iflag = 0
    for i in range(nx):
        if ireq[i] == 2:  #Option 228
            iflag = 2
        if ireq[i] == 3:
            iflag = 3     #Option 229

    #Creating new .fil file if necessary
    if fwhm < 0.0:  #Creating .fil file in case it is necessary
        if iflag==0:
            copyfile(runname+'.fil', path+'/'+runname+'.fil')
            os.chdir(path)
        if iflag==2:  #Option 228
            os.chdir(path)
            iils1 = np.where(ireq == 2)
            iils = iils1[0]
            par1 = xn1[iils[0]]
            par2 = xn1[iils[1]]
            par3 = xn1[iils[2]]
            par4 = xn1[iils[3]]
            par5 = xn1[iils[4]]
            par6 = xn1[iils[5]]
            par7 = xn1[iils[6]]
            ll = ns.write_fil_acsmir(runname,nconv,vconvnem,par1,par2,par3,par4,par5,par6,par7)
        if iflag==3:  #Option 229
            os.chdir(path)
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
    else:
        os.chdir(path)

    xn1 = np.zeros([mx])
    xn1[0:nx] = xn[0:nx]

    #Creating files for communicating with CIRSdrv_wavePY
    xmap = ns.gsetradl(runname,nconv,vconvnem,fwhm,ispace,iscat,gasgiant,layht,nlayer,laytyp,layint,xlat,xlon,lin,hcorrx,nvar,varident,varparam,nx,xn1,jpre,tsurf,isolocc,ionpeel,jlevlo,jlevhi)

   
    #In solar occultation we do not mind about the surface emissivity
    nem = 2
    vem = np.zeros([nem])
    emissivity = np.zeros([nem])
    vem[0]=-100.0
    vem[1]=1.0e7
    emissivity[0:nem]=1.0

    #Writing .cdr file for communicating with CIRSdrv_wavePY
    npath = ngeom
    ll = write_cdr_nemesis(runname,dist,fwhm,ispace,ilbl,nwave,vwave,npath,nconv,vconv,nem,vem,emissivity,tsurf) 

    #Running CIRSdrv_wavePY
    npath,nwave1,nconv1,wavecalc,waveconv,specret_noconv,specret = CIRSdrv_wavePY(runname)
    os.chdir('../')

    #Removing folder
    shutil.rmtree(path)

    return specret_noconv,specret




###############################################################################################


def nemesisSOfm_gas(ix,runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn):

    """
        FUNCTION NAME : nemesisSOfm_gas()
        
        DESCRIPTION : The main purpose from this function is to calculate a forward model, similar to nemesisSOfm, but it is
                      possible to calculate the contribution from just one gas to the spectrum 

                      Reads some input variables like the state vector, and prepares the required files for a CIRSdrv_waveSO run.
                      It also reads an index, which indicates the gas in the .lls which is to be used (0 indicates all gases). Finally,
                      it calls CIRSdrv_waveSO and waits for the output
        
        INPUTS :
        
            ix :: Gas in the .lls to calculate the forward model (0 indicates all gases)
            runname :: Name of the Nemesis run
            nconv :: Number of convolution wavelengths
            vconv(nconv) :: Wavenumber array (cm-1)
            nwave :: Number of calculation wavelenghts
            vwave(nwave) :: Calculation wavenumber array (cm-1)
            npro :: Number of altitude levels in atmosphere
            height(npro) :: Altitude (km)
            ngeom :: Number of tangent heights
            tanhe(ngeom) :: Tangent height (km)
            fwhm :: Full width at half maximum (cm-1) 
            ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
            ilbl :: Flag indicating whether to use correlated-k (0) or lbl (2)
            xlat :: Latitude of observation
            xlon :: Longitude of observation
            lin :: Unit number of previous retrieval (if any)
            nvar :: Number of variables to retrieve
            varident(nvar,3) :: Variable ID
            varparam(nvar,5) :: Other parameters defining the retrieved variables 
            jpre :: Position of tangent pressure (if retrieved)
            nx :: Number of elements in state vector
            xn(nx) :: State vector 

        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            specret(max(nconv),npath) :: Modelled spectra (convolved)
        
        CALLING SEQUENCE:
        
            specret = nemesisSOfm_gas(ix,runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    import nemesisf as ns
    from shutil import copyfile
    from time import sleep

    #Parameters assumed for this version of Nemesis
    iscat = 0 #No scattering
    gasgiant = False #Assumed that it is Mars
    isolocc = 2 #Transmission (not thermal emission)
    tsurf = 100.0 #arbitrary
    dist = 1.5 #arbitrary in transmission
    ionpeel = 1
    layht = 0.0
    laytyp = 5
    layint = 1
    hcorrx = 0.0
    nlayer = npro-1

    #Making array for specifying element in the state vector with special requirements
    ireq = checkvar_nemesisSO(nx,nvar,npro,varident,varparam)

    #Calculating the number oh paths that need to be computed based on 
    mintanhe = tanhe.min()
    maxtanhe = tanhe.max()
    jlevlo1 = np.where(height == mintanhe)
    jlevhi1 = np.where(height == maxtanhe)

    if (len(jlevlo1) > 1):
        sys.exit('error in jacobian_nemesisSO :: jlevlo has more than one value')
    if (len(jlevhi1) > 1):
        sys.exit('error in jacobian_nemesisSO :: jlevhi has more than one value')

    jlevlo = int(jlevlo1[0]) + 1  #because of the fortran indexing
    jlevhi = int(jlevhi1[0]) + 1


    #Making new arrays so that they have the same dimensions as in the Fortran routines
    mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()
    xn1 = np.zeros([mx])
    xn1[0:nx] = xn[:]
    vconvnem = np.zeros([mconv])
    vconvnem[0:nconv] = vconv[:]
    varident1 = np.zeros([mvar,3],dtype='int')
    varparam1 = np.zeros([mvar,mparam])
    varident1[0:nvar,:] = varident[0:nvar,:]
    nparam = len(varparam[0,:])
    varparam1[0:nvar,0:nparam] = varparam[0:nvar,0:nparam]

    #Creating directory to calculate the run and moving necessary files
    path = str(ix)
    mkdir_p(path)

    #Copying necessary files
    copyfile(runname+'.cia', path+'/'+runname+'.cia')
    copyfile(runname+'.ref', path+'/'+runname+'.ref')
    copyfile(runname+'.spx', path+'/'+runname+'.spx')
    copyfile(runname+'.set', path+'/'+runname+'.set')
    copyfile(runname+'.fla', path+'/'+runname+'.fla')
    copyfile(runname+'.sha', path+'/'+runname+'.sha')
    copyfile(runname+'.xsc', path+'/'+runname+'.xsc')
    copyfile(runname+'.lls', path+'/'+runname+'.lls')
    copyfile('aerosol.ref', path+'/'+'aerosol.ref')
    copyfile('height.lay', path+'/'+'height.lay')


    #Checking if ILS is going to be retrieved
    iflag = 0
    for i in range(nx):
        if ireq[i] == 2:  #Option 228
            iflag = 2
        if ireq[i] == 3:
            iflag = 3     #Option 229

    if fwhm < 0.0:  #Creating .fil file in case it is necessary
        if iflag==0:
            copyfile(runname+'.fil', path+'/'+runname+'.fil')
            os.chdir(path)
        if iflag==2:  #Option 228
            os.chdir(path)
            iils1 = np.where(ireq == 2)
            iils = iils1[0]
            par1 = xn1[iils[0]]
            par2 = xn1[iils[1]]
            par3 = xn1[iils[2]]
            par4 = xn1[iils[3]]
            par5 = xn1[iils[4]]
            par6 = xn1[iils[5]]
            par7 = xn1[iils[6]]
            ll = ns.write_fil_acsmir(runname,nconv,vconvnem,par1,par2,par3,par4,par5,par6,par7)
        if iflag==3:  #Option 229
            os.chdir(path)
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
    else:
        os.chdir(path)


    #Selecting the gas in the .lls file to be used for the forward model
    ngasact,strlta = read_lls_nemesis(runname)
    if ix < 0 or ix > ngasact:
        sys.exit('nemesisSOfm_gas :: error. ix must be > 0 and < ngasact')
    
    if ix > 0:
        ngasact = 1
        f = open(runname+'.lls','w')
        f.write(strlta[ix-1])
        f.close()

    xmap = ns.gsetradl(runname,nconv,vconvnem,fwhm,ispace,iscat,gasgiant,layht,nlayer,laytyp,layint,xlat,xlon,lin,hcorrx,nvar,varident1,varparam1,nx,xn1,jpre,tsurf,isolocc,ionpeel,jlevlo,jlevhi)


    #In solar occultation we do not mind about the surface emissivity
    nem = 2
    vem = np.zeros([nem])
    emissivity = np.zeros([nem])
    vem[0]=-100.0
    vem[1]=1.0e7
    emissivity[0:nem]=1.0

    #Writing .cdr file for communicating with CIRSdrv_wavePY
    npath = ngeom
    ll = write_cdr_nemesis(runname,dist,fwhm,ispace,ilbl,nwave,vwave,npath,nconv,vconv,nem,vem,emissivity,tsurf)

    #Running CIRSdrv_wavePY
    npath,nwave1,nconv1,wavecalc,waveconv,specret_noconv,specret = CIRSdrv_wavePY(runname)
    os.chdir('../')

    #Removing folder
    shutil.rmtree(path)

    return specret


###############################################################################################


def nemesisfm_parallel(ip,runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,ngeom,nav,nconv,vconv,fwhm,ispace,gasgiant,\
              layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,flat,flon,lin,nvar,varident,varparam,\
              nx,nxn,xnx,jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwave,nem,vem,emissivity,nproc,xgeom,xav,xx):

    """
        FUNCTION NAME : nemesisfm_parallel()
        
        DESCRIPTION : The main purpose from this function is to calculate a forward model using CIRSdrv_wavePY.
                      Reads some input variables like the state vector, and prepares the required files for a CIRSdrv_wavePY run.
        
        INPUTS :
       
            ip :: Integer indicating which setup for the forward model must be chosen 
            runname :: Name of the Nemesis run
            iscat :: 0=thermal emission; 1=plane parallel scattering; 2= limb/near-limb scattering; 3= single-scattering calculations
            nmu :: Number of zenith ordinates
            mu(nmu) :: Cos(zenith angles)
            wtmu(nmu) :: Quadrature weights
            isol :: Sunlight on(1)/off(0)
            dist :: Solar distance (AU)
            lowbc :: LOwer boundary condition
            galb :: ground albedo
            nf :: Required number of Fourier components
            ngeom :: Number of geometries
            nav(ngeom) :: Number of spectra in each geometry that are required to reconstruct the field-of-view
            nconv(ngeom) :: Number of convolution wavelengths
            vconv(nconv,ngeom) :: Convolution wavelength array for each geometry
            fwhm :: FWHM of convolved spectrum
            ispace :: 0=cm-1, 1= wavelengths
            gasgiant :: Indicates if planet is a gas giant
            layht :: Altitude of base layer (if non-limb, this is set to -emiss_ang for limb calculations)
            nlayer :: Number of atmospheric layers
            laytyp :: How layers are separated
            layint :: How layer amounts are integrated
            sol_ang :: Solar zenith angle
            emiss_ang :: Thermal emission angle
            azi_ang :: Azimuth angle
            flat(ngeom,nav) :: latitude of observation for each geometry and averaging spectrum
            flon(ngeom,nav) :: Longitude of observation for each geometry and averaging spectrum
            lin :: Integer to indicate role of previous retrieval (if any)
            nvar :: Number of retrieved variables
            varident(nvar,3) :: identity of constituent to retrieved and parameterisation
            varparam(nvar,3) :: Additional parameters constraining profile
            nx :: Number of elements in state vector
            nxn :: Number of state vectors that there are in xnx
            xnx(nx,nxn) :: State vectors 
            jalb :: Position of surface albedo spectrum in xn
            jxsc :: Position of x-section spectrum in xn
            jtan :: Position of tangent height correction in xn
            jpre :: Position of tangent pressure
            tsurf :: Surface temperature
            ilbl :: Indicates if correlated-k (0) or lbl (2) is to be used
            nwave(ngeom) :: Number of calculation wavenumbers
            vwave(ngeom,nav) :: Calculation wavenumbers/wavelengths
            nem :: Number of points in surface emissivity spectrum
            vem(nem) :: Wavenumbers for describing the surface emissivity spectrum
            emissivity(nem) :: Surface emissivity spectrum
            nproc :: Total number of processes that are run in parallel
            xgeom(nproc) :: Array indicating which geometry must be chosen for each process to be run
            xav(nproc) :: Array indicating which averagin spectrum (of NAV) must be chosen for each process to be run
            xx(nproc) :: Array indicating which state vector in xnx to use

        OPTIONAL INPUTS: none

        OUTPUTS :

            spec_conv(nconv) :: Modelled spectrum convolved with PSF

        CALLING SEQUENCE:

            spec_conv = nemesisfm_parallel(ip,runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,ngeom,nav,nconv,vconv,fwhm,ispace,gasgiant,
              layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,flat,flon,lin,nvar,varident,varparam,nx,nxn,xnx,
              jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwave,nem,vem,emissivity,nproc,xgeom,xav,xx)


        MODIFICATION HISTORY : Juan Alday (10/06/2019)

    """

    from shutil import copyfile

    ##############################################################################
    #Creating directory to calculate the forward model and moving necessary files
    ##############################################################################

    npath = 1  #We calculate just one path for the forward model in this version (as required by gsetrad)
    path = str(ip)
    mkdir_p(path)

    #Copying necessary files
    copyfile(runname+'.cia', path+'/'+runname+'.cia')
    copyfile(runname+'.ref', path+'/'+runname+'.ref')
    copyfile(runname+'.fla', path+'/'+runname+'.fla')
    copyfile(runname+'.xsc', path+'/'+runname+'.xsc')
    copyfile('aerosol.ref', path+'/'+'aerosol.ref')

    if laytyp==5:   #Base altitude of atmospheric layers must be read from height.lay
        copyfile('height.lay', path+'/'+'height.lay')

    if ilbl==2:   #LBL
        copyfile(runname+'.lls', path+'/'+runname+'.lls')
        if fwhm > 0.0:
            copyfile(runname+'.sha', path+'/'+runname+'.sha')
    elif ilbl==0:  #Correlated-k
        copyfile(runname+'.kls', path+'/'+runname+'.kls')


    if fwhm < 0.0:
        copyfile(runname+'.fil', path+'/'+runname+'.fil')
        os.chdir(path)
    else:
        os.chdir(path)


    ##############################################################################
    #Calling gsetrad to create the rest of the required files
    ##############################################################################

    nconv1 = nconv[xgeom[ip]]
    vconv1 = np.zeros([nconv1])
    vconv1[0:nconv1] = vconv[0:nconv1,xav[ip]]
    nwave1 = nconv[xgeom[ip]]
    vwave1 = np.zeros([nwave1])
    vwave1[0:nwave1] = vconv[0:nwave1,xav[ip]]
    sol_ang1 = sol_ang[xgeom[ip],xav[ip]]
    emiss_ang1 = emiss_ang[xgeom[ip],xav[ip]]
    azi_ang1 = azi_ang[xgeom[ip],xav[ip]]
    xlat = flat[xgeom[ip],xav[ip]]
    xlon = flon[xgeom[ip],xav[ip]]
    xn1 = np.zeros([nx])
    xn1[0:nx] = xnx[0:nx,xx[ip]]
    
    xmap = gsetrad_nemesis(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,\
             layht,nlayer,laytyp,layint,sol_ang1,emiss_ang1,azi_ang1,xlat,xlon,lin,nvar,varident,varparam,nx,xn1,jalb,jxsc,jtan,jpre,tsurf)


    ##############################################################################
    #Writing .cdr file for CIRSrad_wave to read all variables required
    ##############################################################################

    ll = write_cdr_nemesis(runname,dist,fwhm,ispace,ilbl,nwave1,vwave1,npath,nconv1,vconv1,nem,vem,emissivity,tsurf)

    ##############################################################################
    #Running CIRSdrv_wavePY to compute the forward model
    ##############################################################################

    npath,nwave,nconv,vwave,vconv,spec_noconv,spec_conv = CIRSdrv_wavePY(runname)
    os.chdir('../')

    #Removing folder
    shutil.rmtree(path)

    return spec_conv



###############################################################################################


def nemesisfm(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,nconv,vconv,fwhm,ispace,gasgiant,\
              layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,xlat,xlon,lin,nvar,varident,varparam,\
              nx,xn,jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwave,nem,vem,emissivity):

    """
        FUNCTION NAME : nemesisfm()
        
        DESCRIPTION : The main purpose from this function is to calculate a forward model using CIRSdrv_wavePY.
                      Reads some input variables like the state vector, and prepares the required files for a CIRSdrv_wavePY run.
        
        INPUTS :
        
            runname :: Name of the Nemesis run
            iscat :: 0=thermal emission; 1=plane parallel scattering; 2= limb/near-limb scattering; 3= single-scattering calculations
            nmu :: Number of zenith ordinates
            mu(nmu) :: Cos(zenith angles)
            wtmu(nmu) :: Quadrature weights
            isol :: Sunlight on(1)/off(0)
            dist :: Solar distance (AU)
            lowbc :: LOwer boundary condition
            galb :: ground albedo
            nf :: Required number of Fourier components
            nconv :: Number of convolution wavelengths
            vconv(nconv) :: Convolution wavelength array
            fwhm :: FWHM of convolved spectrum
            ispace :: 0=cm-1, 1= wavelengths
            gasgiant :: Indicates if planet is a gas giant
            layht :: Altitude of base layer (if non-limb, this is set to -emiss_ang for limb calculations)
            nlayer :: Number of atmospheric layers
            laytyp :: How layers are separated
            layint :: How layer amounts are integrated
            sol_ang :: Solar zenith angle
            emiss_ang :: Thermal emission angle
            azi_ang :: Azimuth angle
            xlat :: latitude of observation
            xlon :: Longitude of observation
            lin :: Integer to indicate role of previous retrieval (if any)
            nvar :: Number of retrieved variables
            varident(nvar,3) :: identity of constituent to retrieved and parameterisation
            varparam(nvar,3) :: Additional parameters constraining profile
            nx :: Number of elements in state vector
            xn(nx) :: State vector
            jalb :: Position of surface albedo spectrum in xn
            jxsc :: Position of x-section spectrum in xn
            jtan :: Position of tangent height correction in xn
            jpre :: Position of tangent pressure
            tsurf :: Surface temperature
            ilbl :: Indicates if correlated-k (0) or lbl (2) is to be used
            nwave :: Number of calculation wavenumbers
            vwave :: Calculation wavenumbers/wavelengths
            nem :: Number of points in surface emissivity spectrum
            vem(nem) :: Wavenumbers for describing the surface emissivity spectrum
            emissivity(nem) :: Surface emissivity spectrum

        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            spec_noconv(nwave) :: Modelled spectrum at the calculation wavelengths 
            spec_conv(nconv) :: Modelled spectrum convolved with PSF
 
        CALLING SEQUENCE:
        
            spec_noconv,spec_conv = nemesisfm(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,nconv,vconv,fwhm,ispace,gasgiant,
              layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,xlat,lin,nvar,varident,varparam,nx,xn,
              jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwave,nem,vem,emissivity)       

 
        MODIFICATION HISTORY : Juan Alday (10/06/2019)
        
    """

    from shutil import copyfile

    ##############################################################################
    #Creating directory to calculate the forward model and moving necessary files
    ##############################################################################

    npath = 1  #We calculate just one path for the forward model in this version (as required by gsetrad)
    path = 'tmp'
    mkdir_p(path)

    #Copying necessary files
    copyfile(runname+'.cia', path+'/'+runname+'.cia')
    copyfile(runname+'.ref', path+'/'+runname+'.ref')
    copyfile(runname+'.fla', path+'/'+runname+'.fla')
    copyfile(runname+'.xsc', path+'/'+runname+'.xsc')
    copyfile('aerosol.ref', path+'/'+'aerosol.ref')

    if laytyp==5:   #Base altitude of atmospheric layers must be read from height.lay
        copyfile('height.lay', path+'/'+'height.lay')

    if ilbl==2:   #LBL
        copyfile(runname+'.lls', path+'/'+runname+'.lls')
        if fwhm > 0.0:
            copyfile(runname+'.sha', path+'/'+runname+'.sha')
    elif ilbl==0:  #Correlated-k
        copyfile(runname+'.kls', path+'/'+runname+'.kls')


    if fwhm < 0.0:
        copyfile(runname+'.fil', path+'/'+runname+'.fil')
        os.chdir(path)
    else:
        os.chdir(path)


    ##############################################################################
    #Calling gsetrad to create the rest of the required files
    ##############################################################################

    xmap = gsetrad_nemesis(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,nconv,vconv,fwhm,ispace,gasgiant,\
             layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,xlat,xlon,lin,nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf)

    ##############################################################################
    #Writing .cdr file for CIRSrad_wave to read all variables required
    ##############################################################################
 
    ll = write_cdr_nemesis(runname,dist,fwhm,ispace,ilbl,nwave,vwave,npath,nconv,vconv,nem,vem,emissivity,tsurf)


    ##############################################################################
    #Running CIRSdrv_wavePY to compute the forward model
    ##############################################################################

    npath,nwave,nconv,vwave,vconv,spec_noconv,spec_conv = CIRSdrv_wavePY(runname)
    os.chdir('../')

    #Removing folder
    shutil.rmtree(path)


    return spec_noconv,spec_conv


###############################################################################################

def jacobian_nemesis(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,ngeom,nav,nconv,vconv,fwhm,ispace,gasgiant,\
              layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,flat,flon,wgeom,lin,nvar,varident,varparam,\
              nx,xn,jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwave,nem,vem,emissivity,nmaxproc,MakePlot=False,SavePlot=False):

    """
        FUNCTION NAME : jacobian_nemesis()
        
        DESCRIPTION : The main purpose of this function is to calculate the jacobian for a nemesis run, calling several
                      times nemesisfm() in parallel
        
        INPUTS :
        
            runname :: Name of the Nemesis run
            iscat :: 0=thermal emission; 1=plane parallel scattering; 2= limb/near-limb scattering; 3= single-scattering calculations
            nmu :: Number of zenith ordinates
            mu(nmu) :: Cos(zenith angles)
            wtmu(nmu) :: Quadrature weights
            isol :: Sunlight on(1)/off(0)
            dist :: Solar distance (AU)
            lowbc :: LOwer boundary condition
            galb :: ground albedo
            nf :: Required number of Fourier components
            ngeom :: Number of observing geometries
            nav(ngeom) ::  For each geometry, nav how many individual spectral calculations need to be
                           performed in order to reconstruct the field of fiew (then this are all averaged)
            nconv(ngeom) :: Number of convolution wavelengths per geometry
            vconv(nconv,ngeom) :: Convolution wavelength array for each geometry
            fwhm :: FWHM of convolved spectrum
            ispace :: 0=cm-1, 1= wavelengths
            gasgiant :: Indicates if planet is a gas giant
            layht :: Altitude of base layer (if non-limb, this is set to -emiss_ang for limb calculations)
            nlayer :: Number of atmospheric layers
            laytyp :: How layers are separated
            layint :: How layer amounts are integrated
            sol_ang(ngeom,nav) :: Solar zenith angle
            emiss_ang(ngeom,nav) :: Thermal emission angle
            azi_ang(ngeom,nav) :: Azimuth angle
            flat(ngeom,nav) :: Latitude of observation for each geometry
            flon(ngeom,nav) :: Longitude of observation for each geometry
            wgeom(ngeom,nav) :: Weights for each spectrum to be averaged in each geometry 
            lin :: Integer to indicate role of previous retrieval (if any)
            nvar :: Number of retrieved variables
            varident(nvar,3) :: identity of constituent to retrieved and parameterisation
            varparam(nvar,3) :: Additional parameters constraining profile
            nx :: Number of elements in state vector
            xn(nx) :: State vector
            jalb :: Position of surface albedo spectrum in xn
            jxsc :: Position of x-section spectrum in xn
            jtan :: Position of tangent height correction in xn
            jpre :: Position of tangent pressure
            tsurf :: Surface temperature
            ilbl :: Indicates if correlated-k (0) or lbl (2) is to be used
            nwave(ngeom) :: Number of calculation wavenumbers
            vwave(nwave,ngeom) :: Calculation wavenumbers/wavelengths
            nem :: Number of points in surface emissivity spectrum
            vem(nem) :: Wavenumbers for describing the surface emissivity spectrum
            emissivity(nem) :: Surface emissivity spectrum
            nmaxproc :: Maximum number of processes that can be run in parallel at the same time

        OPTIONAL INPUTS:

	    MakePlot :: If True, it makes a plot for the jacobian matrix and the measurement vector
            SavePlot :: If True, saves the plots into .eps files

        OUTPUTS :

            ny :: Number of elements in measurement vector
            yn(ny) :: Measurement vector (same as specret but resized)
            kk(ny,nx) :: Jacobian matrix

        CALLING SEQUENCE:
        
            ny,yn,kk = jacobian_nemesis(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,galb,nf,ngeom,nav,nconv,vconv,fwhm,ispace,gasgiant,\
                                        layht,nlayer,laytyp,layint,sol_ang,emiss_ang,azi_ang,flat,flon,wgeom,lin,nvar,varident,varparam,\
                                        nx,xn,jalb,jxsc,jtan,jpre,tsurf,ilbl,nwave,vwave,nem,vem,emissivity,nmaxproc)

 
        MODIFICATION HISTORY : Juan Alday (10/06/2019)

    """

    import multiprocessing
    from functools import partial
    from shutil import copyfile

    #################################################################################
    # Making some calculations for the proper storage of all the arrays
    #################################################################################

    #For each geometry (and every averaging spectrum) and for each element of the state vector, we need to run a forward model
    nproc = (nx+1) * np.sum(nav)

    #Constructing state vector after perturbation of each of the elements and storing in matrix
    nxn = nx+1
    xnx = np.zeros([nx,nxn])
    for i in range(nx+1):
        if i==0:   #First element is the normal state vector
            xnx[0:nx,i] = xn[0:nx]
        else:      #Perturbation of each element
            xnx[0:nx,i] = xn[0:nx]
            xnx[i-1,i] = xn[i-1]*1.01
            if xn[i-1]==0.0:
                xnx[i-1,i] = 0.05

    #Creating arrays which tell the forward model which geometry and state vector to use
    xgeom = np.zeros([nproc],dtype='int')
    xav = np.zeros([nproc],dtype='int')
    xx = np.zeros([nproc],dtype='int')

    ix = 0
    for i in range(nx+1):
        xx[ix:ix+np.sum(nav)] = i
        for j in range(ngeom):
            xgeom[ix:ix+nav[j]] = j
            for k in range(nav[j]):
                xav[ix+k] = k
            ix = ix + nav[j]


    #Initializing arrays to store all the forward models
    specret_all = np.zeros([nconv.max(),ngeom,nav.max(),nx+1])


    #################################################################################
    # Calling all the forward models
    #################################################################################
    
    #Calculating the array indices that need to be calculated in each set of runs
    nrun = int(nproc/nmaxproc)     #Number of set of runs 
    nrunf = float(nproc)/float(nmaxproc)
    if nrun==0:
        nrun = 1

    if (nrunf - float(nrun)) > 0.0:
        nrun = nrun + 1

    nprocrun = np.zeros([nrun],dtype='int')  #Number of parallel forward models per run
    ix = 0
    for i in range(nrun):
        nprocrun[i] = nmaxproc
        ix = ix + nprocrun[i]
        if ix > nproc:
            ix = ix - nprocrun[i]
            nprocrun[i] = nproc - ix

    
    #Calculating all the forward models
    ix = 0
    for i in range(nrun):
        iproc = range(nprocrun[i])
        for j in range(nprocrun[i]):
            iproc[j] = iproc[j] + ix
        pool = multiprocessing.Pool(processes=nprocrun[i])

        fmproc=partial(nemesisfm_parallel,runname=runname,iscat=iscat,nmu=nmu,mu=mu,wtmu=wtmu,isol=isol,dist=dist,lowbc=lowbc,galb=galb,nf=nf,ngeom=ngeom,\
                       nav=nav,nconv=nconv,vconv=vconv,fwhm=fwhm,ispace=ispace,gasgiant=gasgiant,layht=layht,nlayer=nlayer,laytyp=laytyp,layint=layint,\
                       sol_ang=sol_ang,emiss_ang=emiss_ang,azi_ang=azi_ang,flat=flat,flon=flon,lin=lin,nvar=nvar,varident=varident,varparam=varparam,\
                       nx=nx,nxn=nxn,xnx=xnx,jalb=jalb,jxsc=jxsc,jtan=jtan,jpre=jpre,tsurf=tsurf,ilbl=ilbl,nwave=nwave,vwave=vwave,nem=nem,vem=vem,\
                       emissivity=emissivity,nproc=nproc,xgeom=xgeom,xav=xav,xx=xx)
        result_list = pool.map(fmproc,iproc)

        for j in range(nprocrun[i]):
            specret2 = result_list[j]
            specret_all[0:nconv[xgeom[ix+j]],xgeom[ix+j],xav[ix+j],xx[ix+j]] = specret2[0:nconv[xgeom[ix+j]],0]

        ix = ix + nprocrun[i]


    #################################################################################
    # Averaging for each geometry
    #################################################################################

    specret_ave = np.zeros([nconv.max(),ngeom,nx+1])
    for i in range(nx+1):
        for j in range(ngeom):
            if nav[j]>1:
                spectot = np.zeros([nconv[j]])
                wgeomtot = 0.0
                for k in range(nav[j]):
                    spectot[0:nconv[j]] = spectot[0:nconv[j]] + wgeom[j,k] * specret_all[0:nconv[j],j,k,i]
                    wgeomtot = wgeomtot + wgeom[j,k]                
                specret_ave[0:nconv[j],j,i] = spectot[0:nconv[j]]
            else:
                specret_ave[0:nconv[j],j,i] = specret_all[0:nconv[j],j,0,i]

    #################################################################################
    # Calculating measurement vector and jacobian matrix
    #################################################################################

    ny = sum(nconv)
    yn = np.zeros([ny])
    kk = np.zeros([ny,nx])
      
    ix = 0 
    for j in range(ngeom):
        yn[ix:ix+nconv[j]] = specret_ave[0:nconv[j],j,0]
        ix = ix + nconv[j]

    yn1 = np.zeros([ny])
    for i in range(nx):
        ix = 0
        for j in range(ngeom):
            yn1[ix:ix+nconv[j]] = specret_ave[0:nconv[j],j,i+1]     
            ix = ix + nconv[j]

        kk[0:ny,i] = (yn1[0:ny]-yn[0:ny])/(xnx[i,i+1]-xn[i])


    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        ax = plt.axes()
        ax.set_xlim(0,ny)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Vector element',**axis_font)
        plt.ylabel('Radiance',**axis_font)
        ix = 0
        for i in range(ngeom):
            im = ax.plot(vconv[0:nconv[i],i],yn[ix:ix+nconv[i]])
            ix = ix + nconv[i]
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fm.eps',dpi=100)


        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        ax = plt.axes()
        ax.tick_params(labelsize=20)
        plt.xlabel('Measurement vector y',**axis_font)
        plt.ylabel('State vector x',**axis_font)
        ax.imshow(np.transpose(kk),cmap='hot')
        ax.set_aspect('auto')
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_kk.eps',dpi=100)

    return ny,yn,kk


###############################################################################################


def jacobian_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn,ifix,nmaxproc,MakePlot=False,SavePlot=False):

    """
        FUNCTION NAME : jacobian_nemesisSO()
        
        DESCRIPTION : This function calculates the Jacobian matrix by calling nx+1 times nemesisSOfm(). This routine is set up so
                      that each forward model is calculated in parallel, increasing the computational speed of the code 
        
        INPUTS :
        
            runname :: Name of the Nemesis run
            nconv :: Number of convolution wavelengths
            vconv(nconv) :: Wavenumber array (cm-1)
            nwave :: Number of calculation wavelengths
            vwave(nwave) :: Calculation wavenumber array (cm-1)
            npro :: Number of altitude levels in atmosphere
            height(npro) :: Altitude (km)
            ngeom :: Number of tangent heights
            tanhe(ngeom) :: Tangent height (km)
            fwhm :: Full width at half maximum (cm-1) 
            ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
            ilbl :: Flag indicating whether to use correlated-k (0) or lbl (2)
            xlat :: Latitude of observation
            xlon :: Longitude of observation
            lin :: Unit number of previous retrieval (if any)
            nvar :: Number of variables to retrieve
            varident(nvar,3) :: Variable ID
            varparam(nvar,5) :: Other parameters defining the retrieved variables 
            jpre :: Position of tangent pressure (if retrieved)
            nx :: Number of elements in state vector
            xn(nx) :: State vector
            ifix(nx) :: Array indicating which elements of the state vector shall be fixed (1)
            nmaxproc :: Maximum number of cores to run at the same time 

        OPTIONAL INPUTS: none
        
        OUTPUTS :

            ny :: Number of elements in measurement vector 
            yn(ny) :: Measurement vector (same as specret but resized)              
            kk(ny,nx) :: Jacobian matrix
 
        CALLING SEQUENCE:
        
            ny,yn,kk = jacobian_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn,ifix,nmaxproc)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """


    import multiprocessing
    from functools import partial
    from shutil import copyfile

    #################################################################################
    # Making some calculations for the proper storage of all the arrays
    #################################################################################

    nproc = nx+1 #Number of times we need to run the forward model 

    #Constructing state vector after perturbation of each of the elements and storing in matrix
    nxn = nx+1
    xnx = np.zeros([nx,nxn])
    for i in range(nx+1):
        if i==0:   #First element is the normal state vector
            xnx[0:nx,i] = xn[0:nx]
        else:      #Perturbation of each element
            xnx[0:nx,i] = xn[0:nx]
            xnx[i-1,i] = xn[i-1]*1.01
            if xn[i-1]==0.0:
                xnx[i-1,i] = 0.05


    #################################################################################
    # Calling the forward model all the times
    #################################################################################

    #Calculating the array indices that need to be calculated in each set of runs
    nrun = int(nproc/nmaxproc)     #Number of set of runs 
    nrunf = float(nproc)/float(nmaxproc)
    if nrun==0:
        nrun = 1

    if (nrunf - float(nrun)) > 0.0:
        nrun = nrun + 1

    nprocrun = np.zeros([nrun],dtype='int')  #Number of parallel forward models per run
    ix = 0
    for i in range(nrun):
        nprocrun[i] = nmaxproc
        ix = ix + nprocrun[i]
        if ix > nproc:
            ix = ix - nprocrun[i]
            nprocrun[i] = nproc - ix

    #Running the forward models
    ix = 0
    ny = ngeom*nconv
    yntot = np.zeros([ny,nproc])   #Array in which to store the forward models
    for i in range(nrun):
        iproc = range(nprocrun[i])
        for j in range(nprocrun[i]):
            iproc[j] = iproc[j] + ix
        pool = multiprocessing.Pool(processes=nprocrun[i])
        fmproc=partial(nemesisSOfm_parallel,runname=runname,nconv=nconv,vconv=vconv,nwave=nwave,vwave=vwave,npro=npro,height=height,ngeom=ngeom,tanhe=tanhe,fwhm=fwhm, \
                                   ispace=ispace,ilbl=ilbl,xlat=xlat,xlon=xlon,lin=lin,nvar=nvar,varident=varident,varparam=varparam,jpre=jpre,nx=nx,nxn=nxn,xnx=xnx)
        result_list = pool.map(fmproc,iproc)
        for j in range(nprocrun[i]):
            specret2_1 = result_list[j]
            specret2 = np.zeros([nconv,ngeom])
            specret2[:,:] = specret2_1[:,0:ngeom]
            yn1 = np.resize(np.transpose(specret2),[ny])   #Modelled measurement vector
            yntot[:,iproc[j]] = yn1[:]

        ix = ix + nprocrun[i]

    #If pressure and temperature from just one gas (ireq=1) is to be retrieved, then we need to calcualte
    #the forward model with just that one gas in the .tls file for the reference atmosphere
    flagpt = 0
    for i in range(nvar):
        if (varident[i,0] == 0 and varident[i,1] == -1):
            flagpt = 1
        if (varident[i,0] == 666 and varident[i,1] == -1):
            flagpt = 1

    if flagpt == 1:
       copyfile(runname+'.lls', runname+'_ref.lls')
       copyfile(runname+'.tls', runname+'.lls')
       spec_noconv,specretpt = nemesisSOfm(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn)
       copyfile(runname+'_ref.lls', runname+'.lls')
       ynpt = np.resize(np.transpose(specretpt),[ny])

       ynrest = np.zeros([ngeom*nconv])
       for i in range(ny):
           ynrest[i] = yntot[i,0]/ynpt[i]

    #Reading which elements of the state vector have special requirements
    ireq = checkvar_nemesisSO(nx,nvar,npro,varident,varparam)

    #Calculating the Jacobian matrix
    kk = np.zeros([ny,nx])
    for i in range(nx):
        xn1 = xn[i] * 1.01
        if xn1==0.0:
            xn1=0.05
        if ifix[i] == 0:
            if ireq[i] == 1:  #If p-T are to be retrieved from just one gas
                kk[:,i] = (yntot[:,i+1]-ynpt[:])/(xn1-xn[i])*ynrest[:]
            else:
                kk[:,i] = (yntot[:,i+1]-yntot[:,0])/(xn1-xn[i])


    #Creating array for the measurement vector
    yn = np.zeros([ny])
    yn[0:ny] = yntot[:,0]

    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        wavemin = vconv.min()
        wavemax = vconv.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('Wavenumber (cm$^{-1}$)',**axis_font)
        plt.ylabel('Transmission',**axis_font)
        colors = plt.cm.jet(np.linspace(0,1,ngeom))
        ix = 0
        for i in range(ngeom):
            im = ax.plot(vconv[0:nconv],yn[ix:ix+nconv],color=colors[i])
            ix = ix + nconv
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_fm.eps',dpi=100)


        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(20,8))
        ax = plt.axes()
        ax.tick_params(labelsize=20)
        plt.xlabel('Measurement vector y',**axis_font)
        plt.ylabel('State vector x',**axis_font)
        ax.imshow(np.transpose(kk),cmap='hot')
        ax.set_aspect('auto')
        plt.grid()
        plt.show()
        if SavePlot == True:
            fig.savefig(runname+'_kk.eps',dpi=100)

    return ny,yn,kk

###############################################################################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

###############################################################################################

def calcnextxn_nemesis(nx,ny,xa,xn,y,yn,dd,aa):

    """
        FUNCTION NAME : calcnextxn_nemesis()
        
        DESCRIPTION : 

            This subroutine performs the optimal estimation retrieval of the
            vector x from a set of measurements y and forward derivative matrix
            kk. The equation solved is (re: p147 of Houghton, Taylor and Rodgers):

                             xn+1 = x0 + dd*(y-yn) - aa*(x0 - xn)       
 
        INPUTS :
        
            nx :: Number of elements in state vector
            ny :: Number of elements in measurement vector
            xa(mx) :: A priori state vector
            xn(mx) :: Retrieved state vector in last iteration
            y(my) :: Measurement vector
            yn(my) :: Fitted value of y at last iteration
            dd(mx,my) :: Gain matrix
            aa(mx,mx) :: Averaging kernels

        OPTIONAL INPUTS: none
        
        OUTPUTS :

            x_out(mx) :: Newly retrieved state vector
 
        CALLING SEQUENCE:
        
            x_out = calcnextxn_nemesis(nx,ny,xa,xn,y,yn,dd,aa)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    m1 = np.zeros([ny,1])
    m1[:,0] = y[0:ny] - yn[0:ny]
    dd1 = np.zeros([nx,ny])
    dd1[0:nx,0:ny] = dd[0:nx,0:ny]

    m2 = np.zeros([nx,1])
    m2[:,0] = xa[0:nx] - xn[0:nx]
    aa1 = np.zeros([nx,nx])
    aa1[0:nx,0:nx] = aa[0:nx,0:nx]

    mp1 = np.matmul(dd1,m1)
    mp2 = np.matmul(aa1,m2)

    x_out = np.zeros([nx])

    for i in range(nx):
        x_out[i] = xa[i] + mp1[i,0] - mp2[i,0]

    return x_out


###############################################################################################

def npvar_nemesis(nvar,npro,varident,varparam):

    """
        FUNCTION NAME : npvar_nemesis()
        
        DESCRIPTION : 

            Function for getting the number of points associated with each variable in Nemesis
 
        INPUTS :
        
            nvar :: Number of variables
            npro :: Number of altitude points in atmospheric profiles
            varident(nvar,3) :: Variable ID
            varparam(nvar,mparam) :: Extra parameters for describing the retrieved variable

        OPTIONAL INPUTS: none
        
        OUTPUTS :

            nxvar(nvar) :: Number of points associated with each variable
 
        CALLING SEQUENCE:
        
            nxvar = npvar_nemesis(nvar,npro,varident,varparam)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

#    if nvar==1:
#        imod = varident[2]

    nxvar = np.zeros([nvar],dtype='int')
    for i in range(nvar):
#        if nvar>1:
        imod = varident[i,2]


        if imod == -1:
            nxvar[i] = npro
        elif imod == 0:
            nxvar[i] = npro
        elif imod == 1:
            nxvar[i] = 2
        elif imod == 2:
            nxvar[i] = 1
        elif imod == 3:
            nxvar[i] = 1
        elif imod == 4:
            nxvar[i] = 3       
        elif imod == 5:
            nxvar[i] = 1
        elif imod == 6:
            nxvar[i] = 2
        elif imod == 7:
            nxvar[i] = 2
        elif imod == 8:
            nxvar[i] = 3
        elif imod == 9:
            nxvar[i] = 3
        elif imod == 10:
            nxvar[i] = 4
        elif imod == 11:
            nxvar[i] = 2
        elif imod == 12:
            nxvar[i] = 3
        elif imod == 13:
            nxvar[i] = 3
        elif imod == 14:
            nxvar[i] = 3
        elif imod == 15:
            nxvar[i] = 3
        elif imod == 16:
            nxvar[i] = 4
        elif imod == 17:
            nxvar[i] = 2
        elif imod == 18:
            nxvar[i] = 2
        elif imod == 19:
            nxvar[i] = 4
        elif imod == 20:
            nxvar[i] = 2
        elif imod == 21:
            nxvar[i] = 2
        elif imod == 22:
            nxvar[i] = 5
        elif imod == 23:
            nxvar[i] = 4
        elif imod == 24:
            nxvar[i] = 3
        elif imod == 25:
            nxvar[i] = int(varparam[i,0])
        elif imod == 26:
            nxvar[i] = 4
        elif imod == 27:
            nxvar[i] = 3
        elif imod == 28:
            nxvar[i] = 1
        elif imod == 228:
            nxvar[i] = 7
        elif imod == 229:
            nxvar[i] = 7
        elif imod == 444:
            nxvar[i] = 1 + 1 + int(varparam[i,0])
        elif imod == 666:
            nxvar[i] = 1
        else:
            sys.exit('error :: varID not included in npvar_nemesis()')  
      
    return nxvar



###############################################################################################

def checkvar_nemesisSO(nx,nvar,npro,varident,varparam):

    """
        FUNCTION NAME : checkvar_nemesisSO()
        
        DESCRIPTION : 

            This function creates an array which indicates if there are any elements in the
            state vector which require special treatment for computing the forward model. These
            models are:

              (1) :: Pressure and temperature retrievals (option [0,-1,0] and [666,-1,666]) indicate
                     that the pressure and temperature retrieval must be calculated from the gases
                     included in the .tls file rather than all the active gases in the .lls file.
                     This option is useful for retrieval p-T profiles from CO2 absorption in martian atmosphere.

              (2) :: Indicates that the instrument line shpae is to be retrieved using model 228. Therefore, 
                     the code requires to create a new .fil file.

              (3) :: Indicates that the instrument line shpae is to be retrieved using model 229. Therefore, 
                     the code requires to create a new .fil file.

 
        INPUTS :
       
            nx :: Number of elements in state vector 
            nvar :: Number of variables
            npro :: Number of altitude points in atmospheric profiles
            varident(nvar,3) :: Variable ID
            varparam(nvar,mparam) :: Extra parameters for describing the retrieved variable

        OPTIONAL INPUTS: none
        
        OUTPUTS :

            ireq(nx) :: Array indicating the requirements, as reference in the description
 
        CALLING SEQUENCE:
        
            ireq = checkvar_nemesisSO(nx,nvar,npro,varident,varparam)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

 
    #Making array for specifying element in the state vector with special requirements
    nxvar = npvar_nemesis(nvar,npro,varident,varparam)
    ireq = np.zeros([nx],dtype='int')
    ixx = 0
    for i in range(nvar):
        if (varident[i,0] == 228):
            ireq[ixx:ixx+nxvar[i]] = 2   #ACS-MIR instrument line shape retrieval 
        elif (varident[i,0] == 229):
            ireq[ixx:ixx+nxvar[i]] = 3   #ACS-MIR instrument line shape retrieval
        elif (varident[i,0] == 0 and varident[i,1] == -1):
            ireq[ixx:ixx+nxvar[i]] = 1   #Temperature retrieval from just 1 gas (read .tls file)
        elif (varident[i,0] == 666 and varident[i,1] == -1):
            ireq[ixx:ixx+nxvar[i]] = 1   #Temperature retrieval from just 1 gas (read .tls file)
        ixx = ixx + nxvar[i]

    return ireq


###############################################################################################

def calc_gain_matrix_nemesis(nx,ny,kk,sa,se):


    """
        FUNCTION NAME : calc_gain_matrix_nemesis()
        
        DESCRIPTION : 

            Calculate gain matrix and averaging kernels. The gain matrix is calculated with
               dd = sx*kk_T*(kk*sx*kk_T + se)^-1    (if nx>=ny)
               dd = ((sx^-1 + kk_T*se^-1*kk)^-1)*kk_T*se^-1  (if ny>nx)

 
        INPUTS :
       
            nx :: Number of elements in state vector 
            ny :: Number of elements in measurement vector
            kk(my,mx) :: Jacobian matrix
            sa(mx,mx) :: A priori covariance matric
            se(my,my) :: Measurement error covariance matrix

        OPTIONAL INPUTS: none
        
        OUTPUTS :

            dd(nx,ny) :: Gain matrix
            aa(nx,nx) :: Averaging kernels
 
        CALLING SEQUENCE:
        
            dd,aa = calc_gain_matrix_nemesis(nx,ny,kk,sa,se)       
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Calculating transpose of kk
    kk1 = np.zeros([ny,nx])
    kk1[:,:] = kk[0:ny,0:nx]
    kt1 = np.transpose(kk1)

    #Calculating gain matrix dd
    if (nx >= ny):
        sa1 = np.zeros([nx,nx])   #Resizing just in case
        sa1[:,:] = sa[0:nx,0:nx]
        se1 = np.zeros([ny,ny])
        se1[:,:] = se[0:ny,0:ny]

        #Multiply sa*kt
        m = np.matmul(sa1,kt1)
    
        #Multiply kk*m so that a = kk*sa*kt
        a = np.matmul(kk1,m)    

        #Add se to a so that b = kk*sa*kt + se
        b = np.add(a,se1)

        #Inverting b so that we calculate c = (kk*sa*kt + se)^(-1)
        c = np.linalg.inv(b)

        #Multiplying sa*kt (m above) to c
        dd = np.matmul(m,c)

    else:
        #Calculate inverse of sa and se
        sa1 = np.zeros([nx,nx])   #Resizing just in case
        sa1[:,:] = sa[0:nx,0:nx]
        se1 = np.zeros([ny,ny])
        se1[:,:] = se[0:ny,0:ny]

        sai = np.linalg.inv(sa1)
        sei = np.linalg.inv(se1) 

        #Calculate kt*sei
        m = np.matmul(kt1,sei)
 
        #Calculate m*kk so that kt*se^(-1)*kk
        a = np.matmul(m,kk1)
 
        #Add sai to a so that b = kt*se^(-1)*kk + sa^(-1)
        b = np.add(sai,a)

        #Invert b so that c = (kt*se^(-1)*kk + sa^(-1))^(-1)
        c = np.linalg.inv(b)

        #Multiplying c by kt*sei (m from before) 
        dd = np.matmul(c,m)

    aa = np.matmul(dd,kk1)

    return dd,aa


###############################################################################################


def forwarderr(runname,ngeom,nconv,vconv,woff):
    
    """
        FUNCTION NAME : forwarderr()

        DESCRIPTION : Subroutine which returns the forward modelling error read in from an external file

        INPUTS :
        
            runname :: Name of the Nemesis run (assumes that file has .err extension)
            ngeom :: Number of observing geometries
            nconv(ngeom) :: Number of convolution wavenumbers for each geometry
            vconv(nconv,ngeom) :: Convolution wavenumbers
            woff :: Any wavenumber offset applied to the spectrum

        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            rerr(nconv,ngeom) :: radiance forward-modelling errors

        CALLING SEQUENCE:
        
            rerr = forwarderr(runname,ngeom,nconv,vconv,woff)

        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Reading runname.err file
    f = open(runname+'.err','r')
    tmp = np.fromfile(f,sep=' ',count=1,dtype='int')
    np = int(tmp[0])
    v1 = np.zeros([np])
    r1 = np.zeros([np])
    for i in range(np):
        tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
        v1[i] = float(tmp[0])
        r1[i] = float(tmp[1])
         
    for igeom in range(ngeom):
        for i in range(nconv[igeom]):
            x = vconv[igeom,i]
            if (x < v1[0] or x > v1[np-1]):
                sys.exit('error in forwarderr() :: wavelength not covered in '+runname+'.err file')

    
    sys.exit('ERROR :: FORWARDERR() FUNCTION NEEDS TO BE FINISHED')


###############################################################################################

def calc_serr_nemesis(nx,ny,sa,se,dd,aa):


    """
        FUNCTION NAME : calc_serr_nemesis()
        
        DESCRIPTION : 

            This subroutine calculates the error covariance matrices after the final iteration has been completed.

            The subroutine calculates the MEASUREMENT covariance matrix according to the 
            equation (re: p130 of Houghton, Taylor and Rodgers) :
               
                                  sm = dd*se*dd_T

            The subroutine calculates the SMOOTHING error covariance matrix according to the equation:
  
                                  sn = (aa-I)*sx*(aa-I)_T  

            The subroutine also calculates the TOTAL error matrix:

                                  st=sn+sm

        INPUTS :
       
            nx :: Number of elements in state vector 
            ny :: Number of elements in measurement vector
            sa(mx,mx) :: A priori covariance matric
            se(my,my) :: Measurement error covariance matrix
            dd(mx,my) :: Gain matrix
            aa(mx,mx) :: Averaging kernels

        OPTIONAL INPUTS: none
        
        OUTPUTS :

            sm(nx,nx) :: Final measurement covariance matrix
            sn(nx,nx) :: Final smoothing error covariance matrix
            st(nx,nx) :: Final full covariance matrix
 
        CALLING SEQUENCE:
        
            sm,sn,st = calc_gain_matrix_nemesis(nx,ny,sa,se,dd,aa)       
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """

    #Resizing matrices just in case
    dd1 = np.zeros([nx,ny])
    dd1[:,:] = dd[0:nx,0:ny]
    aa1 = np.zeros([nx,nx])
    aa1[:,:] = aa[0:nx,0:nx]
    sa1 = np.zeros([nx,nx])
    sa1[:,:] = sa[0:nx,0:nx]
    se1 = np.zeros([ny,ny])
    se1[:,:] = se[0:ny,0:ny]

    dt1 = np.transpose(dd1)
  
    #Multiplying dd*se
    a = np.matmul(dd1,se1)

    #Multiplying a*dt so that dd*se*dt
    sm = np.matmul(a,dt1)

    #Calculate aa-ii where I is a diagonal matrix
    b = np.zeros([nx,nx])
    for i in range(nx):
        for j in range(nx):
            b[i,j] = aa[i,j]
        b[i,i] = b[i,i] - 1.0
    bt = np.transpose(b)

    #Multiply b*sa so that (aa-I)*sa
    c = np.matmul(b,sa1)
  
    #Multiply c*bt so tthat (aa-I)*sx*(aa-I)_T  
    sn = np.matmul(c,bt)

    #Add sn and sm and get total retrieved error
    st = np.add(sn,sm)

    return sm,sn,st


###############################################################################################


def calc_phiret_nemesis(ny,y,yn,se,nx,xn,xa,sa):


    """
        FUNCTION NAME : calc_phiret_nemesis()
        
        DESCRIPTION : 

            Calculate the retrieval cost function which combines departure from a priori and closeness to spectrum.
 
        INPUTS :
      
            ny :: Number of elements in measurement vector
            y(my) :: Measurement vector
            yn(my) :: Modelled measurement vector
            se(my,my) :: Measurement error covariance matrix
            nx :: Number of elements in state vector 
            xn(mx) :: State vector
            xa(mx) :: A priori state vector
            sa(mx,mx) :: A priori covariance matrix       
 
        OPTIONAL INPUTS: none
        
        OUTPUTS :

            chisq :: Closeness of fit to measurement vector
            phi :: Total cost function
 
        CALLING SEQUENCE:
        
            chisq,phi = calc_phiret_nemesis(ny,y,yn,se,nx,xn,xa,sa)       
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Calculating yn-y
    b = np.zeros([ny,1])
    b[:,0] = yn[0:ny] - y[0:ny]
    bt = np.transpose(b)

    #Calculating inverse of sa and se
    sa1 = np.zeros([nx,nx])   #Resizing just in case
    sa1[:,:] = sa[0:nx,0:nx]
    se1 = np.zeros([ny,ny])
    se1[:,:] = se[0:ny,0:ny]

    sai = np.linalg.inv(sa1)
    sei = np.linalg.inv(se1)

    #Multiplying se^(-1)*b
    a = np.matmul(sei,b)
 
    #Multiplying bt*a so that (yn-y)^T * se^(-1) * (yn-y)
    c = np.matmul(bt,a)

    phi1 = c[0,0]
    chisq = phi1

    #Calculating xn-xa
    d = np.zeros([nx,1])
    d[:,0] = xn[0:nx] - xa[0:nx]
    dt = np.transpose(d)
   
    #Multiply sa^(-1)*d 
    e = np.matmul(sai,d)

    #Multiply dt*e so that (xn-xa)^T * sa^(-1) * (xn-xa)
    f = np.matmul(dt,e)

    phi2 = f[0,0]
   
    print('calc_phiret_nemesis: phi1,phi2 = '+str(phi1)+','+str(phi2)+')')
    phi = phi1+phi2

    return chisq,phi



###############################################################################################


def assess_nemesis(nx,ny,kk,sa,se):


    """
        FUNCTION NAME : assess_nemesis()
        
        DESCRIPTION : 

            This subroutine assesses the retrieval matrices to see
            whether an exact retrieval may be expected.

            One formulation of the gain matrix is dd = sx*kk_T*(kk*sx*kk_T + se)^-1

            If the retrieval is exact, the se will be very small. Since se is
            diagonal all we need do is compare to  the diagonal elements of
 
        INPUTS :
      
            nx :: Number of elements in state vector 
            ny :: Number of elements in measurement vector
            kk(my,mx) :: Jacobian matrix
            sa(mx,mx) :: A priori covariance matric
            se(my,my) :: Measurement error covariance matrix
 
        OPTIONAL INPUTS: none
        
        OUTPUTS : none

        CALLING SEQUENCE:
        
            dummy = assess_nemesis(nx,ny,kk,sa,se)       
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Calculating transpose of kk
    kk1 = np.zeros([ny,nx])
    kk1[:,:] = kk[0:ny,0:nx]
    kt1 = np.transpose(kk1)

    #Multiply sa*kt
    sa1 = np.zeros([nx,nx])   #Resizing just in case
    sa1[:,:] = sa[0:nx,0:nx]
    se1 = np.zeros([ny,ny])
    se1[:,:] = se[0:ny,0:ny]
    m = np.matmul(sa1,kt1)

    #Multiply kk*m so that a = kk*sa*kt
    a = np.matmul(kk1,m)

    #Add se to a
    b = np.add(a,se1)

    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    for i in range(ny):
            sum1 = sum1 + b[i,i]
            sum2 = sum2 + se1[i,i]
            sum3 = sum3 + b[i,i]/se1[i,i]

    sum1 = sum1/ny
    sum2 = sum2/ny
    sum3 = sum3/ny
  
    print('Assess:')
    print('Average of diagonal elements of Kk*Sx*Kt : '+str(sum1))
    print('Average of diagonal elements of Se : '+str(sum2))
    print('Ratio = '+str(sum1/sum2))
    print('Average of Kk*Sx*Kt/Se element ratio : '+str(sum3))
    if sum3 > 10.0:
        print('******************* ASSESS WARNING *****************')
        print('Insufficient constraint. Solution likely to be exact')
        print('****************************************************')

    dummy = 0
    return dummy


###############################################################################################


def write_fil_acsmir_v2(runname,nconv,vconv,par1,par2,par3,par4,par5,par6,par7):


    """
        FUNCTION NAME : write_fil_acsmir_v2()
        
        DESCRIPTION : 

            Write the .fil file using the approach of a double gaussian for ACS-MIR spectra
 
        INPUTS :
      
            runname :: Name of the Nemesis run
            nconv :: Number of convolution wavelengths
            wave(nconv) :: Convolution wavelengths
            par1 :: Wavenumber offset of main at lowest wavenumber
            par2 :: Wavenumber offset of main at wavenumber in the middle
            par3 :: Wavenumber offset of main at highest wavenumber 
            par4 :: Offset of the second gaussian with respect to the first one (assumed spectrally constant)
            par5 :: FWHM of the main gaussian at lowest wavenumber (assumed to be constat in wavelength units)
            par6 :: Relative amplitude of the second gaussian with respect to the gaussian at lowest wavenumber
            par7 :: Relative amplitude of the second gaussian with respect to the gaussian at highest wavenumber (linear var)

        OPTIONAL INPUTS: none
        
        OUTPUTS : 

            runname.file file

        CALLING SEQUENCE:
        
            ll = write_fil_acsmir_v2(runname,nconv,vconv,par1,par2,par3,par4,par5,par6,par7)       
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    print('WRITE_FIL_ACSMIR :: Calculating the ILS')
    print('Parameters :: ')
    print('Par. 1 = '+str(par1))
    print('Par. 2 = '+str(par2))
    print('Par. 3 = '+str(par3))
    print('Par. 4 = '+str(par4))
    print('Par. 5 = '+str(par5))
    print('Par. 6 = '+str(par6))
    print('Par. 7 = '+str(par7))

    #Opening .fil file
    f = open(runname+'.fil','w')
    f.write("%i \n" %  (nconv))
    

    #Calculating the parameters for each spectral point
    vconv1 = np.zeros([nconv])    
    vconv1[:] = vconv[0:nconv]
    ng = 2

    # 1. Wavenumber offset of the two gaussians
    #    We divide it in two sections with linear polynomials     
    iconvmid = int(nconv/2.)
    wavemax = vconv1[nconv-1]
    wavemin = vconv1[0]
    wavemid = vconv1[iconvmid]
    offgrad1 = (par2 - par1)/(wavemid-wavemin)
    offgrad2 = (par2 - par3)/(wavemid-wavemax)
    offset = np.zeros([nconv,ng])
    for i in range(iconvmid):
        offset[i,0] = (vconv1[i] - wavemin) * offgrad1 + par1
        offset[i,1] = offset[i,0] + par4
    for i in range(nconv-iconvmid):
        offset[i+iconvmid,0] = (vconv1[i+iconvmid] - wavemax) * offgrad2 + par3
        offset[i+iconvmid,1] = offset[i+iconvmid,0] + par4

    # 2. FWHM for the two gaussians
    fwhm = np.zeros([nconv,ng])
    fwhml = par5 / wavemin**2.0
    for i in range(nconv):
        fwhm[i,0] = fwhml * (vconv[i])**2.
        fwhm[i,1] = fwhm[i,0]

    # 3. Amplitde of the second gaussian with respect to the main one
    amp = np.zeros([nconv,ng])
    ampgrad = (par7 - par6)/(wavemax-wavemin)
    for i in range(nconv):
        amp[i,0] = 1.0
        amp[i,1] = (vconv[i] - wavemin) * ampgrad + par6

    
    #Running for each spectral point
    for i in range(nconv):
        f.write("%10.7f\n" % vconv[i])


        #determining the lowest and highest wavenumbers to calculate
        xlim = 0.0
        xdist = 5.0 
        for j in range(ng):
            xcen = offset[i,j]
            xmin = abs(xcen - xdist*fwhm[i,j]/2.)
            if xmin > xlim:
                xlim = xmin
            xmax = abs(xcen + xdist*fwhm[i,j]/2.)
            if xmax > xlim:
                xlim = xmax


        #determining the wavenumber spacing we need to sample properly the gaussians
        xsamp = 7.0   #number of points we require to sample one HWHM 
        xhwhm = 10000.0
        for j in range(ng):
            xhwhmx = fwhm[i,j]/2. 
            if xhwhmx < xhwhm:
                xhwhm = xhwhmx
        deltawave = xhwhm/xsamp
        np1 = 2.0 * xlim / deltawave
        npx = int(np1) + 1

        #Calculating the ILS in this spectral point
        iamp = np.zeros([ng])
        imean = np.zeros([ng])
        ifwhm = np.zeros([ng])
        fun = np.zeros([npx])
        xwave = np.linspace(vconv[i]-deltawave*(npx-1)/2.,vconv[i]+deltawave*(npx-1)/2.,npx)        
        for j in range(ng):
            iamp[j] = amp[i,j]
            imean[j] = offset[i,j] + vconv[i]
            ifwhm[j] = fwhm[i,j]

        fun = ngauss(npx,xwave,ng,iamp,imean,ifwhm)  

        #Writing output
        f.write("%i \n" %  (npx))
        for j in range(npx):
            f.write("%10.10f %10.10e\n" % (xwave[j], fun[j]) )


    f.close()

    dummy = 0.0
    return dummy


###############################################################################################


def ngauss(npx,x,ng,iamp,imean,ifwhm,MakePlot=False):


    """
        FUNCTION NAME : ngauss()
        
        DESCRIPTION : 

            Create a function which is the sum of multiple gaussians
 
        INPUTS :
      
            npx :: Number of points in x-array
            x(npx) :: Array specifying the points at which the function must be calculated
            ng :: Number of gaussians
            iamp(ng) :: Amplitude of each of the gaussians
            imean(ng) :: Center x-point of the gaussians
            ifwhm(ng) :: FWHM of the gaussians

        OPTIONAL INPUTS: none
        
        OUTPUTS : 

            fun(npx) :: Function at each x-point

        CALLING SEQUENCE:
        
            fun = ngauss(npx,x,ng,iamp,imean,ifwhm)
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    fun  = np.zeros([npx])
    isigma = ifwhm/(2.0*np.sqrt(2.*np.log(2.)))
    for i in range(npx):
        for j in range(ng):
            fun[i] = fun[i] + iamp[j] * np.exp( -(x[i]-imean[j])**2.0/(2.0*isigma[j]**2.)  )


    #Make plot if keyword is specified
    if MakePlot == True:
        axis_font = {'fontname':'Arial', 'size':'20'}
        cm = plt.cm.get_cmap('RdYlBu')
        fig = plt.figure(figsize=(15,8))
        wavemin = x.min()
        wavemax = x.max()
        ax = plt.axes()
        ax.set_xlim(wavemin,wavemax)
        ax.tick_params(labelsize=20)
        ax.ticklabel_format(useOffset=False)
        plt.xlabel('x',**axis_font)
        plt.ylabel('f(x)',**axis_font)
        im = ax.plot(x,fun)
        plt.grid()
        plt.show()    
    
    return fun


###############################################################################################


def write_cov_nemesis(runname,npro,nvar,varident,varparam,nx,ny,sa,sm,sn,st,se,aa,dd,kk):


    """
        FUNCTION NAME : write_cov_nemesis()
        
        DESCRIPTION : 

            Reads the input parameters and writes the .cov file with the standard Nemesis format
 
        INPUTS :
      
            runname :: Name of the Nemesis run 
            npro :: Number of points in atmospheric profiles
            nvar :: Number of retrieved variables
            varident(nvar,3) :: Variable ID
            varparam(nvar,mparam) :: Extra parameters for describing the retrieved variable 
            nx :: Number of elements in state vector 
            ny :: Number of elements in measurement vector
            sa(mx,mx) :: A priori covariance matric
            sm(nx,nx) :: Final measurement covariance matrix
            sn(nx,nx) :: Final smoothing error covariance matrix
            st(nx,nx) :: Final full covariance matrix
            se(my,my) :: Measurement error covariance matrix
            aa(mx,mx) :: Averaging kernels
            dd(mx,my) :: Gain matrix
            kk(my,mx) :: Jacobian matrix

        OPTIONAL INPUTS: none
        
        OUTPUTS : 

            Nemesis .cov file

        CALLING SEQUENCE:
        
            dummy = write_cov_nemesis(runname,npro,nvar,varident,varparam,nx,ny,sa,sm,sn,st,se,aa,dd,kk)
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """


    #Open file
    f = open(runname+'.cov','w')

    f.write("%i %i\n" % (npro,nvar))

    for i in range(nvar):
        f.write("%i \t %i \t %i\n" % (varident[i,0],varident[i,1],varident[i,2]))
        f.write("%10.8e \t %10.8e \t %10.8e \t %10.8e \t %10.8e\n" % (varparam[i,0],varparam[i,1],varparam[i,2],varparam[i,3],varparam[i,4]))

    f.write("%i %i\n" % (nx,ny))

    for i in range(nx):
        for j in range(nx):
            f.write("%10.8e\n" % (sa[i,j]))
        for j in range(nx):
            f.write("%10.8e\n" % (sm[i,j]))
        for j in range(nx):
            f.write("%10.8e\n" % (sn[i,j]))
        for j in range(nx):
            f.write("%10.8e\n" % (st[i,j]))

    for i in range(nx):
        for j in range(nx):
            f.write("%10.8e\n" % (aa[i,j]))

    for i in range(nx):
        for j in range(ny):
            f.write("%10.8e\n" % (dd[i,j]))

    for i in range(ny):
        for j in range(nx):
            f.write("%10.8e\n" % (kk[i,j]))

    for i in range(ny):
        f.write("%10.8e\n" % (se[i,i]))

    f.close() 

    dummy = 1
    return dummy

###############################################################################################

def write_mre_nemesis(runname,iform,ispace,xlat,xlon,npro,nvar,varident,varparam,nx,ny,y,yn, \
                      se,xa,sa,xn,st,ngeom,nconv,vconv,gasgiant,jpre,jrad,jlogg,iscat,lin,nspec):


    """
        FUNCTION NAME : write_mre_nemesis()
        
        DESCRIPTION : 

            Reads the inputs and writes the .mre file as given by the standard Nemesis format
 
        INPUTS :
      
            runname :: Name of the Nemesis run 
            iform :: Output format 
                      0 = radiance
                      1 = F_plan/F_star 
                      2 = 100*A_plan/A_star
                      3 = planet spectral flux
                      4 = Transmission*solar_flux
                      5 = Transmission (solar occultation)
            ispace :: (0) Wavenumber (cm-1) (1) Wavelength (um)
            xlat :: Latitude
            xlon :: Longitude
            npro :: Number of points in atmospheric profiles
            nvar :: Number of variables to retrieve
            varident(nvar,3) :: Identity of profiles and parametrisation
            varparam(nvar,nparam) :: Extra parameters as required
            nx :: Number of elements in state vector
            ny :: Number of elements in measurement vector
            y(ny) :: Measured spectrum
            yn(ny) :: Best calculated spectrum
            se(ny,ny) :: Measurement covariance matrix
            xa(nx) :: A priori state vector
            sa(nx,nx) :: A priori covariance matrix
            xn(nx) :: Retrieved state vector
            st(nx,nx) :: Retrieved error covariance matrix
            ngeom :: Number of observation geometries
            nconv(ngeom) :: Number of convolution wavenumbers at each geometry
            vconv(nconv,ngeom) :: Convolution wavenumbers
            gasgiant :: Gas giant flag
            jpre :: Indicates if pressure retrieval performed
            jrad :: Indicates if radiusretrieval performed
            jlogg :: Indicates if surface gravity retrieval performed
            iscat :: Flag to indicate scattering calc.
            lin :: Previous retrieval flag

        OPTIONAL INPUTS: none
        
        OUTPUTS : 

            Nemesis .cov file

        CALLING SEQUENCE:
        
            dummy = write_mre_nemesis(runname,iform,ispace,xlat,xlon,npro,nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,st,ngeom,nconv,vconv,gasgiant,jpre,jrad,jlogg,iscat,lin,nspec)
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """
   
    import nemesisf as ns
  
    #Reading the array sizes that the fortran arrays require
    mx,my,mconv,mwave,maxpat,maxlay,maxgas,maxsec,mgeom,mvar,mparam,maxmu = check_arraysize_nemesis()    

    #Opening file
    f = open(runname+'.mre','w')
    
    str1 = '! Total number of retrievals'
    f.write("\t" + str(nspec)+ "\t" + str1 + "\n")

    for ispec in range(nspec):
 
        #Writing first lines
        ispec1 = ispec + 1
        str2 = '! ispec,ngeom,ny,nx,ny'
        f.write("\t %i %i %i %i %i \t %s \n" % (ispec,ngeom,ny,nx,ny,str2)) 
        str3 = 'Latitude, Longitude'
        f.write("\t %5.7f \t %5.7f \t %s \n" % (xlat,xlon,str3)) 
    
        if ispace==0: #Wavenumber space
            if iform==0:
                str4='Radiances expressed as nW cm-2 sr-1 cm'       
                xfac=1.0e9
            elif iform==1:
                str4='F_plan/F_star Ratio of planet'
                xfac=1.0
            elif iform==3:
                str4='Spectral Radiation of planet: W (cm-1)-1'
                xfac=1.0e18
            elif iform==2:
                str4='Transit depth: 100*Planet_area/Stellar_area'
                xfac=1.0
            elif iform==4:
                str4='Solar flux: W cm-2 (cm-1)-1'
                xfac=1.0
            elif iform==5:
                str4='Transmission'
                xfac=1.0
            else:
                print('write_mre_nemesis :: ERROR! iform not defined. Default=0')
                str4='Radiances expressed as nW cm-2 sr-1 cm' 
                xfac=1.0e9
        elif ispace==1: #Wavelength space
            if iform==0:
                str4='Radiances expressed as uW cm-2 sr-1 um-1'
                xfac=1.0e6
            elif iform==1:
                str4='F_plan/F_star Ratio of planet'
                xfac=1.0
            elif iform==3:
                str4='Spectral Radiation of planet: W um-1'
                xfac=1.0e18
            elif iform==2:
                str4='Transit depth: 100*Planet_area/Stellar_area'
                xfac=1.0
            elif iform==4:
                str4='Solar flux: W cm-2 um-1'
                xfac=1.0
            elif iform==5:
                str4='Transmission'
                xfac=1.0
            else:
                print('write_mre_nemesis :: ERROR! iform not defined. Default=0')
                str4='Radiances expressed as uW cm-2 sr-1 um-1'
                xfac=1.0e6

        f.write(str4+"\n")

        #Writing spectra
        l = ['i','lambda','R_meas','error','%err','R_fit','%Diff']
        f.write("\t %s %s %s %s %s %s %s \n" % (l[0],l[1],l[2],l[3],l[4],l[5],l[6]))
        ioff = 0
        for igeom in range(ngeom):
            for iconv in range(nconv[igeom]):
                i = ioff+iconv
                err1 = np.sqrt(se[i,i])
                if y[i] != 0.0:
                    xerr1 = abs(100.0*err1/y[i])
                    relerr = abs(100.0*(y[i]-yn[i])/y[i])
                else:
                    xerr1=-1.0
                    relerr1=-1.0

                if iform==0:
                    strspec = "\t %4i %14.8f %15.8e %15.8e %7.2f %15.8f %9.5f \n"
                elif iform==1:
                    strspec = "\t %4i %10.4f %15.8e %15.8e %7.2f %15.8f %9.5f \n"
                elif iform==2:
                    strspec = "\t %4i %9.4f %12.6e %12.6e %6.2f %12.6f %6.2f \n"
                elif iform==3:
                    strspec = "\t %4i %10.4f %15.8e %15.8e %7.2f %15.8f %9.5f \n"
                else:
                    strspec = "\t %4i %14.8f %15.8e %15.8e %7.2f %15.8f %9.5f \n"

                f.write(strspec % (i+1,vconv[iconv,igeom],y[i]*xfac,err1*xfac,xerr1,yn[i]*xfac,relerr))
                
            ioff = ioff + nconv[igeom]      

        #Writing a priori and retrieved state vectors
        str1 = '' 
        f.write(str1+"\n")
        nvar1 = nvar
        f.write('nvar=    '+str(nvar1)+"\n")
        #Reading number of points associated with each variable
        nxvar = npvar_nemesis(nvar,npro,varident,varparam)

        #Reading the variables are in logscale or not
        logvar = logflag_nemesis(nvar,npro,varident,varparam)
        nxtemp = 0
        for ivar in range(nvar):
            f.write('Variable '+str(ivar+1)+"\n")
            f.write("\t %i \t %i \t %i\n" % (varident[ivar,0],varident[ivar,1],varident[ivar,2]))
            f.write("%10.8e \t %10.8e \t %10.8e \t %10.8e \t %10.8e\n" % (varparam[ivar,0],varparam[ivar,1],varparam[ivar,2],varparam[ivar,3],varparam[ivar,4]))
            
            l = ['i','ix','xa','sa_err','xn','xn_err']
            f.write("\t %s %s %s %s %s %s\n" % (l[0],l[1],l[2],l[3],l[4],l[5]))
            for ip in range(nxvar[ivar]):
                ix = nxtemp + ip 
                xa1 = xa[ix]
                ea1 = np.sqrt(abs(sa[ix,ix]))
                xn1 = xn[ix]
                en1 = np.sqrt(abs(st[ix,ix]))
                if logvar[ix]==1:
                    xa1 = np.exp(xa1)
                    ea1 = xa1*ea1
                    xn1 = np.exp(xn1)
                    en1 = xn1*en1

                strx = "\t %4i %4i %12.5e %12.5e %12.5e %12.5e \n"
                f.write(strx % (ip+1,ix+1,xa1,ea1,xn1,en1))
            nxtemp = nxtemp + nxvar[ivar]    
            
    f.close()


###############################################################################################


def write_heightlay_nemesis(nlayer,heightlay):


    """
        FUNCTION NAME : write_heightlay_nemesis()
        
        DESCRIPTION : 

            Writes the height.lay file with the input required by Nemesis. This file specifies the
            base altitude of each layer in the atmosphere
 
        INPUTS :
      
            nlayer :: Number of layers in atmosphere
            heightlay(nlayer) :: Base altitude of each layer (km)

        OPTIONAL INPUTS: none
        
        OUTPUTS : 

            Nemesis height.lay file

        CALLING SEQUENCE:
        
            dummy = write_heightlay_nemesis(nlayer,heightlay)
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    f = open('height.lay','w')
    header = 'Nemesis simulation - base altitude of atmospheric layers'
    f.write(header+"\n")
    f.write('\t %i \n' % (nlayer))
    for i in range(nlayer):
        f.write('\t %7.2f \n' % (heightlay[i]))

    f.close()


###############################################################################################


def calcgascn_nemesisSO(runname,nconv,vconv,nwave,vwave,npro,height,ngeom,tanhe,fwhm,ispace,ilbl,xlat,xlon,lin,nvar,varident,varparam,jpre,nx,xn,WRITE_GCN=False):


    """
        FUNCTION NAME : calcgascn_nemesisSO()
        
        DESCRIPTION : 

            After a retrieval has been performed, it calculates the contribution from each of the
            active gases, calling several times nemesisSOfm.
 
        INPUTS :

            runname :: Name of the Nemesis run
            nconv :: Number of convolution wavelengths
            vconv(nconv) :: Wavenumber array (cm-1)
            nwave :: Number of calculation wavelengths
            vwave(nwave) :: Calculation wavenumber array (cm-1)
            npro :: Number of altitude levels in atmosphere
            height(npro) :: Altitude (km)
            ngeom :: Number of tangent heights
            tanhe(ngeom) :: Tangent height (km)
            fwhm :: Full width at half maximum (cm-1) 
            ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
            ilbl :: Flag indicating whether to use correlated-k (0) or lbl (2)
            xlat :: Latitude of observation
            xlon :: Longitude of observation
            lin :: Unit number of previous retrieval (if any)
            nvar :: Number of variables to retrieve
            varident(nvar,3) :: Variable ID
            varparam(nvar,5) :: Other parameters defining the retrieved variables 
            jpre :: Position of tangent pressure (if retrieved)
            nx :: Number of elements in state vector
            xn(nx) :: State vector      

        OPTIONAL INPUTS: 

	    WRITE_GCN :: If True, the function writes a file with the same output as the nemesisSO .gcn file,
                         which can be read using read_gcn_nemesisSO()
        
        OUTPUTS : 

            specret(nconv,ngeom) :: Spectra with the contribution from ALL active gases
            specretgas(nconv,ngeom,ngas) :: Spectra with the contribution from each for the active gases

        CALLING SEQUENCE:
        
            specret,specretgas = calcgascn_nemesisSO(runname,nconv,vconv,npro,height,ngeom,tanhe,fwhm,ispace,xlat,lin,nvar,varident,varparam,jpre,nx,xn)
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """


    import multiprocessing
    from functools import partial
    from shutil import copyfile


    ngasact,strlta = read_lls_nemesis(runname)
    nproc = ngasact + 1  #Number of times we need to run the forward model

    iproc = range(nproc)
    pool = multiprocessing.Pool(processes=nproc)
    fmproc=partial(nemesisSOfm_gas,runname=runname,nconv=nconv,vconv=vconv,nwave=nwave,vwave=vwave,npro=npro,height=height,ngeom=ngeom,tanhe=tanhe,fwhm=fwhm, \
                               ispace=ispace,ilbl=ilbl,xlat=xlat,xlon=xlon,lin=lin,nvar=nvar,varident=varident,varparam=varparam,jpre=jpre,nx=nx,xn=xn)
    result_list = pool.map(fmproc,iproc)

    specret = np.zeros([nconv,ngeom])
    specretgas = np.zeros([nconv,ngeom,ngasact])
    for j in range(nproc):
        specret2 = result_list[j]
        if j == 0:
            specret[:,:] = specret2[:,0:ngeom]    
        else:
            specretgas[:,:,j-1] = specret2[:,0:ngeom]


    if WRITE_GCN == True:
        inst_fwhm,xlat,xlon,ngeom,nav,wgeom,nconv1,flat,flon,tanhe,wave,meas,errmeas = read_spx_nemesisl(runname)
        gasIDs = np.zeros(ngasact,dtype='int')
        isoIDs = np.zeros(ngasact,dtype='int')
        for i in range(ngasact):
            nwave,vmin,delv,npress,ntemp,gasID,isoID,presslevels,templevels = read_ltahead_nemesis(strlta[i])
            gasIDs[i] = gasID
            isoIDs[i] = isoID

        f = open(runname+'.gcn','w')
        f.write('\t %i \t %i \t %i \n' % (nconv,ngeom,ngasact))
        for i in range(ngasact):
            f.write('\t %i \t %i \n' % (gasIDs[i],isoIDs[i]))
        for i in range(ngeom):
            f.write('\t %10.6f \n' % (tanhe[i]))
            for j in range(nconv):
                str1 = str('{0:7.6f}'.format(vconv[j]))+'\t'+str('{0:7.6f}'.format(meas[j,i,0]))+'\t'+str('{0:7.6f}'.format(errmeas[j,i,0]))+'\t'+str('{0:7.6f}'.format(specret[j,i]))
                for k in range(ngasact):
                    str1 = str1+'\t'+str('{0:7.6e}'.format(specretgas[j,i,k]))
                f.write(str1+'\n')
        f.close()

    return specret,specretgas



###############################################################################################


def read_apr_nemesis(runname,npro):


    """
        FUNCTION NAME : read_apr_nemesis()
        
        DESCRIPTION : 

            Reads the .apr file, which contains information about the variables and 
            parametrisations that are to be retrieved, as well as their a priori values
 
            N.B. In this code, the apriori and retrieved vectors x are usually 
            converted to logs, all except for temperature and fractional scale
            heights
            This is done to reduce instabilities when different parts of the
            vectors and matrices hold vastly different sized properties. e.g. 
            cloud x-section and base height.

        INPUTS :
      
            runname :: Name of the Nemesis run
            npro :: Number of elements in atmospheric profiles

        OPTIONAL INPUTS: none
        
        OUTPUTS : 

            nvar :: Number of variables
            varident(nvar,3) :: Variable ID as presented in the Nemesis manual
            varparam(nvar,nparam) :: Additional parameters constraining the profile
            jsurf :: Position of surface temperature element (if included)
            jalb :: Position of start of surface albedo spectrum (if included)
            jxsc :: Position of start of x-section spectrum (if included)
            jtan :: Position of tangent altitude correction (if included)
            jpre :: Position of ref. tangent  pressure (if included)
            jrad :: Position of radius of planet (if included)
            jlogg :: Position of surface log_10(g) of planet (if included)
            nx :: Number of elements in state vector
            xa(nx) :: A priori state vector
            sa(nx,nx) :: A priori covariance matrix
            lx(nx) :: Log flag. 0 if real number 1 if log number 
            csx(nvar) :: Ratio volume of shell/total volume of particle for Maltmieser 
                         coated sphere model. For homogeneous sphere model, csx(ivar)=-1


        CALLING SEQUENCE:
        
            nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,nx,xa,sa,lx,csx = read_apr_nemesis(runname)
 
        MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Open file
    f = open(runname+'.apr','r')

    #Reading header
    s = f.readline().split()

    #Reading first line
    tmp = np.fromfile(f,sep=' ',count=1,dtype='int')
    nvar = int(tmp[0])

    #Initialise some variables
    jsurf = -1
    jalb = -1
    jxsc = -1
    jtan = -1
    jpre = -1
    jrad = -1
    jlogg = -1
    sxminfac = 0.001    
    mparam = 200        #Giving big sizes but they will be re-sized
    mx = 1000 
    varident = np.zeros([nvar,3],dtype='int')
    varparam = np.zeros([nvar,mparam])
    lx = np.zeros([mx],dtype='int')
    x0 = np.zeros([mx])
    sx = np.zeros([mx,mx])

    #Reading data
    ix = 0
    for i in range(nvar):
        tmp = np.fromfile(f,sep=' ',count=3,dtype='int')
        varident[i,:] = tmp[:]

        #Starting different cases
        if varident[i,2] <= 100:    #Parameter must be an atmospheric one
            
            if varident[i,2] == 0: 
#           ********* continuous profile ************************
                s = f.readline().split()
                f1 = open(s[0],'r')
                tmp = np.fromfile(f1,sep=' ',count=2,dtype='float')
                nlevel = int(tmp[0])
                if nlevel != npro:
                    sys.exit('profiles must be listed on same grid as .prf')
                clen = float(tmp[1])
                pref = np.zeros([nlevel])
                ref = np.zeros([nlevel])
                eref = np.zeros([nlevel]) 
                for j in range(nlevel):
                    tmp = np.fromfile(f1,sep=' ',count=3,dtype='float')
                    pref[j] = float(tmp[0])
                    ref[j] = float(tmp[1])
                    eref[j] = float(tmp[2])
                f1.close()
          
                if varident[i,0] == 0:  # *** temperature, leave alone ****
                    x0[ix:ix+nlevel] = ref[:]
                    for j in range(nlevel):
                        sx[ix+j,ix+j] = eref[j]**2.
                else:                   #**** vmr, cloud, para-H2 , fcloud, take logs ***
                    for j in range(nlevel):
                        lx[ix+j] = 1
                        x0[ix+j] = np.log(ref[j])
                        sx[ix+j,ix+j] = ( eref[j]/ref[j]  )**2. 
 
                #Calculating correlation between levels in continuous profile
                for j in range(nlevel):
                    for k in range(nlevel):
                        if pref[j] < 0.0:
                            sys.exit('Error in read_apr_nemesis().  A priori file must be on pressure grid')

                        delp = np.log(pref[k])-np.log(pref[j])
                        arg = abs(delp/clen)
                        xfac = np.exp(-arg)
                        if xfac >= sxminfac:
                             sx[ix+j,ix+k] = np.sqrt(sx[ix+j,ix+j]*sx[ix+k,ix+k])*xfac
                             sx[ix+k,ix+j] = sx[ix+j,ix+k]

                ix = ix + nlevel

            elif varident[i,2] == -1:
#           * continuous cloud, but cloud retrieved as particles/cm3 rather than 
#           * particles per gram to decouple it from pressure.
#           ********* continuous particles/cm3 profile ************************
                if varident[i,0] >= 0:
                    sys.exit('error in read_apr_nemesis :: model -1 type is only for use with aerosols')

                s = f.readline().split()
                f1 = open(s[0],'r')
                tmp = np.fromfile(f1,sep=' ',count=2,dtype='float')
                nlevel = int(tmp[0])
                if nlevel != npro:
                    sys.exit('profiles must be listed on same grid as .prf')
                clen = float(tmp[1])
                pref = np.zeros([nlevel])
                ref = np.zeros([nlevel])
                eref = np.zeros([nlevel])
                for j in range(nlevel):
                    tmp = np.fromfile(f1,sep=' ',count=3,dtype='float')
                    pref[j] = float(tmp[0])
                    ref[j] = float(tmp[1])
                    eref[j] = float(tmp[2])
 
                    lx[ix+j] = 1
                    x0[ix+j] = np.log(ref[j])
                    sx[ix+j,ix+j] = ( eref[j]/ref[j]  )**2.

                f1.close()

                #Calculating correlation between levels in continuous profile
                for j in range(nlevel):
                    for k in range(nlevel):
                        if pref[j] < 0.0:
                            sys.exit('Error in read_apr_nemesis().  A priori file must be on pressure grid')

                        delp = np.log(pref[k])-np.log(pref[j])
                        arg = abs(delp/clen)
                        xfac = np.exp(-arg)
                        if xfac >= sxminfac:
                             sx[ix+j,ix+k] = np.sqrt(sx[ix+j,ix+j]*sx[ix+k,ix+k])*xfac
                             sx[ix+k,ix+j] = sx[ix+j,ix+k]

                ix = ix + nlevel


            elif varident[i,2] == 1:
#           ******** profile held as deep amount, fsh and knee pressure ** 
#           Read in xdeep,fsh,pknee
                tmp = np.fromfile(f,sep=' ',count=1,dtype='float') 
                pknee = float(tmp[0])
                tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
                xdeep = float(tmp[0])
                edeep = float(tmp[1])
                tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
                xfsh = float(tmp[0])
                efsh = float(tmp[1])

                varparam[i,0] = pknee
                
                if varident[i,0] == 0:  #Temperature, leave alone
                    x0[ix] = xdeep
                    sx[ix,ix] = edeep**2.
                else:
                    x0[ix] = np.log(xdeep)
                    sx[ix,ix] = ( edeep/xdeep )**2.
                    lx[ix] = 1

                ix = ix + 1
               
                if xfsh > 0.0:
                    x0[ix] = np.log(xfsh)
                    lx[ix] = 1
                    sx[ix,ix] = ( efsh/xfsh  )**2.
                else:
                    sys.exit('Error in read_apr_nemesis().  xfsh must be > 0')
 
                ix = ix + 1



            elif varident[i,2] == 2:
#           **** Simple scaling factor of reference profile *******
#           Read in scaling factor
 
                tmp = np.fromfile(f,sep=' ',count=2,dtype='float')              
                x0[ix] = float(tmp[0])
                sx[ix,ix] = (float(tmp[1]))**2.

                ix = ix + 1

            elif varident[i,2] == 3:
#           **** Exponential scaling factor of reference profile *******
#           Read in scaling factor
 
                tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
                xfac = float(tmp[0])
                err = float(tmp[1])

                if xfac > 0.0:
                    x0[ix] = np.log(xfac)
                    lx[ix] = 1
                    sx[ix,ix] = ( err/xfac ) **2.
                else:
                    sys.exit('Error in read_apr_nemesis().  xfac must be > 0')

                ix = ix + 1


            elif varident[i,2] == 4:
#           ******** profile held as deep amount, fsh and VARIABLE knee press
#           Read in xdeep,fsh,pknee
                tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
                pknee = float(tmp[0])
                eknee = float(tmp[1])
                tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
                xdeep = float(tmp[0])
                edeep = float(tmp[1])
                tmp = np.fromfile(f,sep=' ',count=2,dtype='float')
                xfsh = float(tmp[0])
                efsh = float(tmp[1]) 


                if varident[i,0] == 0:  #Temperature, leave alone
                    x0[ix] = xdeep
                    sx[ix,ix] = edeep**2.
                else:
                    x0[ix] = np.log(xdeep)
                    sx[ix,ix] = ( edeep/xdeep )**2.
                    lx[ix] = 1
                ix = ix + 1

                if xfsh > 0.0:
                    x0[ix] = np.log(xfsh)
                    lx[ix] = 1
                    sx[ix,ix] = ( efsh/xfsh  )**2.
                else:
                    sys.exit('Error in read_apr_nemesis().  xfsh must be > 0')
                ix = ix + 1

                x0[ix] = np.log(pknee)
                lx[ix] = 1
                sx[ix,ix] = (eknee/pknee)**2
                ix = ix + 1

            else:
                sys.exit('error in read_apr_nemesis() :: Variable ID not included in this function')





