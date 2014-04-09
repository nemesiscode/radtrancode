      subroutine readapriori(opfile,lin,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,nx,x0,sx,lx)
C     $Id:
C     ****************************************************************
C     Subroutine to read in apriori vector and covariance matrix
C
C     Input variables
C	opfile		character*100	root name of run.
C	lin		integer		Indicates if previous retrieval
C					to be used as a priori for new
C					retrieval
C	lpre		integer		Previous retrieval unit number
C	xlat		real		Required latitude
C	npro		integer		Number if vertical levels in .prf
C
C     Output variables
C	nvar		integer		Number of variable profiles
C					(including T, vmr and cloud)
C	varident(mvar,3)integer		identity of constituent to 
C  					retrieved and how it is represented
C					First and second column contains
C					identity. Third column contains:
C					0 read in new profile and error
C					1 read in deep, fsh, knee
C					2 scale profile in .ref file.
C       varparam(mvar,mparam) integer   Additional parameters constraining
C					 profile.
C	jsurf		integer		Position of surface temperature 
C					element (if included)
C 	jalb		integer		Position of start of surface
C					albedo spectrum
C	jtan		integer		Position of tangent altitude
C					correction
C	jpre		integer		Position of ref. tangent  pressure
C       jrad            integer         Position of radius of planet
C       jlogg           integer         Position of surface log_10(g) of planet
C	nx 		integer 	number of elements in state vector
C	x0(mx)		real		a priori vector
C	sx(mx,mx)	real		a priori covariance matrix
C	lx(mx)		integer		Log flag. 0 if real number
C						  1 if log number 
C
C     N.B. In this code, the apriori and retrieved vectors x are usually 
C     converted to logs, all except for temperature and fractional scale
C     heights
C     This is done to reduce instabilities when different parts of the
C     vectors and matrices hold vastly different sized properties. e.g. 
C     cloud x-section and base height.
C 
C     Original:	Pat Irwin		10/10/03
C     31/10/03 PGJI Revised to limit the smallest off-diagonal element in 
C                   sx to be greater than sx(i,i)*sx(j,j)*SXMINFAC. This
C                   was done because when inverting sx before, 
C                   instabilities crept in when sx contained very 
C		    small, but not quite zero off-diagonal elements!
C
C     ****************************************************************

      implicit none

      integer i,j,nx,ix,jx,npro,jsurf,np,jalb,jtan,jpre,jrad,maxlat,k
      integer jlogg,nmode,nwave,max_mode, max_wave
      parameter (max_mode = 10)
      parameter (max_wave = 1000)
      parameter(maxlat=100)

C     ****************************************************************
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/planrad.f'
      include 'arraylen.f'
C     ****************************************************************

      real xsec(max_mode,max_wave,2),wave(max_wave)
      real x0(mx),sx(mx,mx),err,err1,ref1,pref(maxpro)
      real eref(maxlat,maxpro),reflat(maxlat),htan
      real delp,xfac,pknee,eknee,edeep,xdeep,xlat,xlatx,xlonx,pshld
      real efsh,xfsh,varparam(mvar,mparam),flat,hknee,pre
      real ref(maxlat,maxpro),clen,SXMINFAC,arg,valb,alb
      real xknee,xrh,erh,xcdeep,ecdeep,radius,Grav
      parameter (Grav=6.672E-11)
C     SXMINFAC is minimum off-diagonal factor allowed in the
C     a priori covariance matrix
      parameter (SXMINFAC = 0.001)
      real varparamx(mvar,mparam),xnx(mx),sxx(mx,mx),xdiff
      integer varident(mvar,3),ivar,nvar,nlevel,lin,jsurfx
      integer jalbx,jtanx,jprex,jradx,jlat,ilat,nlat,lx(mx)
      integer nprox,nvarx,varidentx(mvar,3),lpre,ioffx,ivarx
      integer npx,ioff,icond,npvar,jloggx,iplanet
      character*100 opfile,buffer,ipfile,runname,rifile
      integer nxx 
      real xwid,ewid,y,y0,lambda0,vi(mx)
      real r0,er0,dr,edr,vm,nm,nimag,delv
      real xldeep,eldeep,xlhigh,elhigh,v1


C     Initialise a priori parameters
      do i=1,mx
       x0(i)=0.0
       lx(i)=0
       do j=1,mx
        sx(j,i)=0.
       end do
      end do
      jsurf = -1
      jalb = -1
      jtan = -1
      jpre = -1
      jrad = -1
      jlogg = -1
      runname=opfile

C     Pass jrad and jlogg to planrad common block
      jradf=jrad
      jloggf=jlogg

      call file(opfile,opfile,'apr')
      open(27,file=opfile,status='old')

      read(27,1)buffer
1     format(a)
    
      read(27,*)nvar			! Number of variable profiles
      print*,'Reading in ',nvar,' variables'
      if(nvar.gt.mvar) then
       print*,'Readapriori: NVAR can not be greater than MVAR'
       print*,nvar,mvar
       stop
      endif

C     varident(ivar,1-2) holds the identity of the profile. The default is gas
C     ID and ISO, e.g. 6 1.
C     if varident(ivar,1) = 0, then the profile holds temperature
C     if varident(ivar,1) < 0, then the profile holds cloud opacity with
C     icont = abs(varident(ivar,1))

C     varident(ivar,3) holds how the profile is to be read in
C     0 read in new profile and error as a function of pressure
C     1 read in deep, fsh, knee
C     2 scale profile in .ref file.
C     3 scale profile in .ref file by the exponent of a number.


      nx = 0 
      do 10 ivar=1,nvar
         read(27,*)(varident(ivar,j),j=1,3)
         if(varident(ivar,1).le.100)then
C          parameter must be an atmospheric one.

           if(varident(ivar,3).eq.0)then
C          ********* continuous profile ************************
             read(27,1)ipfile
             print*,'reading variable ',ivar,' from ',ipfile
             open(28,file=ipfile,status='old')
             read(28,*)nlevel,clen
             if(nlevel.ne.npro)then
              print*,'profiles must be listed on same grid as .prf'
              stop
             endif
             do 20 i=1,nlevel
              read(28,*)pref(i),ref(1,i),eref(1,i)
              ix = i+nx
C             For vmrs and cloud density always hold the log. 
C             Avoids instabilities arising from greatly different profiles 
C             such as T and vmrs
              if(varident(ivar,1).eq.0)then
C              *** temperature, leave alone ****
               x0(ix) = ref(1,i)
               err = eref(1,i)
              else
C              **** vmr, cloud, para-H2 , fcloud, take logs ***
               if(ref(1,i).gt.0.0) then
                 x0(ix) = alog(ref(1,i)) 
                 lx(ix)=1
               else 
                 print*,'Error in readapriori.f. Cant take log of zero'
                 print*,i,ref(i,i)
                 stop
               endif
               err = eref(1,i)/ref(1,i)
              endif
              sx(ix,ix)=err**2

20           continue
             close(28)
             do i = 1,nlevel
              ix = nx + i         
              do j = i+1,nlevel
               jx = nx + j
               if(pref(i).lt.0.0) then
                 print*,'Error in readapriori.f. A priori file '
                 print*,'must be on pressure grid '
                 stop
               endif            
c               print*,'DIAG1',pref(i),pref(j),alog(pref(j))
               delp = log(pref(j))-log(pref(i))
c               print*,'168: CLEN:',clen
c               print*,'169: DELP:',delp
               
               arg = abs(delp/clen)
               xfac = exp(-arg)
C               xfac = exp(-arg*arg)
               if(xfac.ge.SXMINFAC)then  
                sx(ix,jx) = sqrt(sx(ix,ix)*sx(jx,jx))*xfac
                sx(jx,ix) = sx(ix,jx)
               endif
              enddo
             enddo

             nx = nx + nlevel

           elseif (varident(ivar,3).eq.1)then
C            ******** profile held as deep amount, fsh and knee pressure ** 
C            Read in xdeep,fsh,pknee
             read(27,*)pknee
             read(27,*)xdeep,edeep
             read(27,*)xfsh,efsh
             varparam(ivar,1) = pknee
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)=xdeep
              err = edeep
             else
C             *** vmr, para-H2 or cloud, take logs *********
              if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
              else
                print*,'Error in readapriori. xdeep must be > 0.0'
                stop
              endif
              err = edeep/xdeep
             endif
             sx(ix,ix)=err**2
             ix = nx+2
             if(xfsh.gt.0.0)then
               x0(ix) = alog(xfsh)
	       lx(ix)=1
             else
               print*,'Error in readapriori - xfsh must be > 0'
               stop
             endif
             sx(ix,ix) = (efsh/xfsh)**2

             nx = nx+2

           elseif(varident(ivar,3).eq.2)then
C            **** Simple scaling factor of reference profile *******
C            Read in scaling factor
             ix = nx+1
             read(27,*)x0(ix),err
             sx(ix,ix) = err**2

             nx = nx+1

           elseif(varident(ivar,3).eq.3)then
C            **** Exponential scaling factor of reference profile *******
C            Read in scaling factor
             ix = nx+1
             read(27,*)xfac,err
             if(xfac.gt.0.0)then
               x0(ix)=alog(xfac)
               lx(ix)=1
             else
               print*,'Error in readpriori - xfac must be > 0'
               stop
             endif
             err = err/xfac
             sx(ix,ix) = err**2

             nx = nx+1

           elseif (varident(ivar,3).eq.4)then
C            ******** profile held as deep amount, fsh and VARIABLE knee press
C            Read in xdeep,fsh,pknee
             read(27,*)pknee,eknee
             read(27,*)xdeep,edeep
             read(27,*)xfsh,efsh
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)=xdeep
              err = edeep
             else
C             *** vmr, fcloud, para-H2 or cloud, take logs *********
              if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
              else
                print*,'Error in readapriori. xdeep must be > 0.0'
                stop
              endif
              err = edeep/xdeep
             endif
             sx(ix,ix)=err**2
             ix = nx+2
             if(xfsh.gt.0.0)then
               x0(ix) = alog(xfsh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xfsh must be > 0'
               stop
             endif
             sx(ix,ix) = (efsh/xfsh)**2
             ix = nx+3
             x0(ix) = alog(pknee)
             lx(ix)=1

             sx(ix,ix) = (eknee/pknee)**2

             nx = nx+3

           elseif(varident(ivar,3).eq.5)then
C           *********** a priori profile is continuous and varies with latitude
             read(27,1)ipfile
             print*,'reading variable ',ivar,' from ',ipfile
             open(28,file=ipfile,status='old')
             read(28,*)nlevel,nlat,clen
             if(nlevel.ne.npro)then
              print*,'profiles must be listed on same grid as .prf'
              stop
             endif
             do 26 ilat=1,nlat
              read(28,*)reflat(ilat)
              do 25 i=1,nlevel
               read(28,*)pref(i),ref(ilat,i),eref(ilat,i)
25            continue
26           continue
             
             jlat=0
             do 37 ilat=1,nlat-1
              if(xlat.ge.reflat(ilat).and.xlat.lt.reflat(ilat+1))then
               flat = (xlat-reflat(ilat))/(reflat(ilat+1)-reflat(ilat))
               jlat = ilat 
              endif
37           continue
           
             if(jlat.lt.1)then
              print*,'Error in readapriori'
              print*,'Apriori file does not cover latitude range'
              print*,nlat
              print*,(reflat(i),i=1,nlat)
              stop
             endif
             
             do 27 i=1,nlevel
              ix = i+nx
C             For vmrs and cloud density always hold the log. 
C             Avoids instabilities arising from greatly different profiles 
C             such as T and vmrs
              ref1 = (1.0-flat)*ref(jlat,i)+flat*ref(jlat+1,i)
              err1 = (1.0-flat)*eref(jlat,i)+flat*eref(jlat+1,i)              
              if(varident(ivar,1).eq.0)then
C              *** temperature, leave alone ****
               x0(ix) = ref1
               err = err1
              else
C              **** vmr, cloud, para-H2 , fcloud, take logs ***
               if(ref1.gt.0)then
                  x0(ix) = alog(ref1)
                  lx(ix)=1
               else
                  print*,'Error in readapriori - ref1 must be > 0'
                  stop
               endif
               err = err1/ref1
              endif
              sx(ix,ix)=err**2

27           continue
             close(28)

             do i = 1,nlevel
              ix = nx + i         
              do j = i+1,nlevel
               jx = nx + j         
               delp = alog(pref(j))-alog(pref(i))
               arg = abs(delp/clen)
               xfac = exp(-arg)
C               xfac = exp(-arg*arg)
               if(xfac.ge.SXMINFAC)then  
                sx(ix,jx) = sqrt(sx(ix,ix)*sx(jx,jx))*xfac
                sx(jx,ix) = sx(ix,jx)
               endif
              enddo
             enddo

             nx = nx + nlevel

           elseif(varident(ivar,3).eq.6)then
C           Venus-type cloud profile. Needs a reference height,
C           integrated optical depth above that height and a cloud
C           scale height (km)

C            Read in xdeep,fsh,pknee
             read(27,*)hknee
             read(27,*)xdeep,edeep
             read(27,*)xfsh,efsh
             varparam(ivar,1) = hknee
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)=xdeep
              err = edeep
             else
C             *** vmr, para-H2 or cloud, take logs *********
              if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
              else
                print*,'Error in readapriori. xdeep must be > 0.0'
                stop
              endif
              err = edeep/xdeep
             endif
             sx(ix,ix)=err**2
             ix = nx+2
             if(xfsh.gt.0.0)then
               x0(ix) = alog(xfsh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xfsh must be > 0'
               stop
             endif
             sx(ix,ix) = (efsh/xfsh)**2

             nx = nx+2

           elseif (varident(ivar,3).eq.7)then
C            ******** profile held as a value at a certain pressure and 
C            ******** fractional scale height 
C            Read in xknee,fsh,pknee
             read(27,*)pknee
             read(27,*)xknee,eknee
             read(27,*)xfsh,efsh
             varparam(ivar,1) = pknee
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)=xknee
              err = eknee
             else
C             *** vmr, fcloud, para-H2 or cloud, take logs *********
              if(xknee.gt.0)then
                x0(ix)=alog(xknee)
                lx(ix)=1
              else
                print*,'Error in readapriori - xknee must be > 0'
                stop
              endif
              err = eknee/xknee
             endif
             sx(ix,ix)=err**2
             ix = nx+2
             if(xfsh.gt.0.0)then
               x0(ix) = alog(xfsh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xfsh must be > 0'
               stop
             endif
             sx(ix,ix) = (efsh/xfsh)**2

             nx = nx+2

           elseif (varident(ivar,3).eq.8)then
C            ******** profile held as value at a VARIABLE knee pressure
C            ******** plus a fractional scale height. Below the knee
C            ******** pressure the profile is set to zero - a simple
C            ******** cloud in other words!
C            Read in xdeep,fsh,pknee
             read(27,*)pknee,eknee
             read(27,*)xdeep,edeep
             read(27,*)xfsh,efsh
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)=xdeep
              err = edeep
             else
C             *** vmr, fcloud, para-H2 or cloud, take logs *********
              if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
              else
                print*,'Error in readapriori. xdeep must be > 0.0'
                stop
              endif
              err = edeep/xdeep
             endif
             sx(ix,ix)=err**2
             ix = nx+2
             if(xfsh.gt.0.0)then
               x0(ix) = alog(xfsh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xfsh must be > 0'
               stop
             endif
             sx(ix,ix) = (efsh/xfsh)**2
             ix = nx+3
             if(pknee.gt.0)then
                x0(ix)=alog(pknee)
                lx(ix)=1
             else
                print*,'Error in readapriori - pknee must be > 0'
                stop
             endif
             sx(ix,ix) = (eknee/pknee)**2

             nx = nx+3

           elseif (varident(ivar,3).eq.9)then
C            ******** cloud profile held as total optical depth plus
C            ******** base height and fractional scale height. Below the knee
C            ******** pressure the profile is set to zero - a simple
C            ******** cloud in other words!
C            Read in xdeep,fsh,pknee
             read(27,*)hknee,eknee
             read(27,*)xdeep,edeep
             read(27,*)xfsh,efsh
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)=xdeep
              err = edeep
             else
C             *** vmr, fcloud, para-H2 or cloud, take logs *********
              if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
              else
                print*,'Error in readapriori. xdeep must be > 0.0'
                stop
              endif
              err = edeep/xdeep
             endif
             sx(ix,ix)=err**2
             ix = nx+2
             if(xfsh.gt.0.0)then
               x0(ix) = alog(xfsh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xfsh must be > 0'
               stop
             endif
             sx(ix,ix) = (efsh/xfsh)**2
             ix = nx+3
             x0(ix) = hknee
             sx(ix,ix) = eknee**2

             nx = nx+3

           elseif (varident(ivar,3).eq.10)then
C            Variable condensible gas and associated cloud
             read(27,*)xdeep,edeep
             read(27,*)xrh,erh
             read(27,*)xcdeep,ecdeep
             read(27,*)xfsh,efsh
             read(27,*)varparam(ivar,1) 

             ix = nx+1
             if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
             else
               print*,'Error in readapriori. xdeep must be > 0.0'
               stop
             endif
             err = edeep/xdeep
             sx(ix,ix)=err**2

             ix = nx+2
             if(xrh.gt.0.0)then
               x0(ix) = alog(xrh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xrh must be > 0'
               stop
             endif
             sx(ix,ix) = (erh/xrh)**2

             ix = nx+3
             if(xcdeep.gt.0.0)then
                x0(ix)=alog(xcdeep)
                lx(ix)=1
             else
                print*,'Error in readapriori. xcdeep must be > 0.0'
                stop
             endif
             err = ecdeep/xcdeep
             sx(ix,ix)=err**2

             ix = nx+4
             if(xfsh.gt.0.0)then
               x0(ix) = alog(xfsh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xfsh must be > 0'
               stop
             endif
             sx(ix,ix) = (efsh/xfsh)**2

             nx=nx+4


           elseif (varident(ivar,3).eq.11)then
C            Variable deep abundance and relative humidity
             read(27,*)xdeep,edeep
             read(27,*)xrh,erh
C            read flag to say if you want the RH to apply at all levels (0)
C                                  or only above the condensation level (1)
             read(27,*)icond

             varparam(ivar,1)=icond

             ix = nx+1
             if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
             else
               print*,'Error in readapriori. xdeep must be > 0.0'
               stop
             endif
             err = edeep/xdeep
             sx(ix,ix)=err**2

             ix = nx+2
             if(xrh.gt.0.0)then
               x0(ix) = alog(xrh)
               lx(ix)=1
             else
               print*,'Error in readapriori - xrh must be > 0'
               stop
             endif
             sx(ix,ix) = (erh/xrh)**2

             nx=nx+2

           elseif (varident(ivar,3).eq.12)then
C            ** profile held as peak amount, pressure and FWHM (Gaussian) **
C            Read in xdeep, pknee, xwid

             read(27,*)xdeep,edeep
             read(27,*)pknee,eknee
             read(27,*)xwid,ewid

             ix = nx+1
             if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
             else
               print*,'Error in readapriori. xdeep must be > 0.0'
               stop
             endif

             err = edeep/xdeep
             sx(ix,ix)=err**2

             ix = nx+2
             x0(ix) = alog(pknee)
             lx(ix)=1

             sx(ix,ix) = (eknee/pknee)**2

             ix = nx+3
             x0(ix) = alog(xwid)
             lx(ix)=1

             sx(ix,ix) = (ewid/xwid)**2

             nx = nx+3


           elseif (varident(ivar,3).eq.13)then
C            ** profile held as peak amount, pressure and FWHM (Lorentzian) **
C            Read in xdeep, pknee, xwid

             read(27,*)xdeep,edeep
             read(27,*)pknee,eknee
             read(27,*)xwid,ewid

             ix = nx+1
             if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
             else
               print*,'Error in readapriori. xdeep must be > 0.0'
               stop
             endif

             err = edeep/xdeep
             sx(ix,ix)=err**2

             ix = nx+2
             x0(ix) = alog(pknee)
             lx(ix)=1

             sx(ix,ix) = (eknee/pknee)**2

             ix = nx+3
             x0(ix) = alog(xwid)
             lx(ix)=1

             sx(ix,ix) = (ewid/xwid)**2

             nx = nx+3


           elseif (varident(ivar,3).eq.14)then
C            ** profile held as peak amount, pressure and FWHM (Gaussian) **
C            Read in xdeep, pknee, xwid

             read(27,*)xdeep,edeep
             read(27,*)hknee,eknee
             read(27,*)xwid,ewid


             ix = nx+1
             if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
             else
               print*,'Error in readapriori. xdeep must be > 0.0'
               stop
             endif

             err = edeep/xdeep
             sx(ix,ix)=err**2

             ix = nx+2
             x0(ix) = hknee
             sx(ix,ix) = eknee**2

             ix = nx+3
             x0(ix) = alog(xwid)
             lx(ix)=1

             sx(ix,ix) = (ewid/xwid)**2

             nx = nx+3


           elseif (varident(ivar,3).eq.15)then
C            ** profile held as peak amount, pressure and FWHM (Lorentzian) **
C            Read in xdeep, pknee, xwid

             read(27,*)xdeep,edeep
             read(27,*)hknee,eknee
             read(27,*)xwid,ewid

             ix = nx+1
             if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
             else
               print*,'Error in readapriori. xdeep must be > 0.0'
               stop
             endif

             err = edeep/xdeep
             sx(ix,ix)=err**2


             ix = nx+2
             x0(ix) = hknee
             sx(ix,ix) = eknee**2

             ix = nx+3
             x0(ix) = alog(xwid)
             lx(ix)=1

             sx(ix,ix) = (ewid/xwid)**2

             nx = nx+3

           elseif (varident(ivar,3).eq.16)then
C            ** profile held as amount at certain pressure with specified
C            lapse rate above and below
C            Read in xdeep, pknee, xwid

             read(27,*)xdeep,edeep
             read(27,*)pknee,eknee
             read(27,*)xldeep,eldeep
             read(27,*)xlhigh,elhigh

             ix = nx+1

             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)=xdeep
              err = edeep
             else
C             *** vmr, para-H2 or cloud, take logs *********
              if(xdeep.gt.0.0)then
                x0(ix)=alog(xdeep)
                lx(ix)=1
              else
                print*,'Error in readapriori. xdeep must be > 0.0'
                stop
              endif
              err = edeep/xdeep
             endif

             sx(ix,ix)=err**2


             ix = nx+2
             x0(ix) = alog(pknee)
             sx(ix,ix) = (eknee/pknee)**2
             lx(ix)=1

             ix = nx+3
             x0(ix) = alog(xldeep)
             lx(ix)=1

             sx(ix,ix) = (eldeep/xldeep)**2

             ix = nx+4
             x0(ix) = alog(xlhigh)
             lx(ix)=1

             sx(ix,ix) = (elhigh/xlhigh)**2

             nx = nx+4

           else         
            print*,'vartype profile parametrisation not recognised'
            stop
           endif
	    
         else
C          Non-atmospheric and other parameters

C          extra line to force varident(ivar,3)=varident(ivar,1)
           varident(ivar,3)=varident(ivar,1)

           if(varident(ivar,1).eq.999)then
C           **** Surface temperature *******
C           Read in surface temperature and error
            ix = nx+1
            read(27,*)x0(ix),err
            sx(ix,ix) = err**2
            jsurf = ix

            nx = nx+1

           elseif(varident(ivar,1).eq.888)then
C           **** Surface albedo spectrum *******
C           Read in number of points
            read(27,*)np
            varparam(ivar,1)=np
            jalb = nx+1
            do i=1,np             
             ix = jalb+i-1
             read(27,*)valb,alb,err
             if(alb.gt.0.0)then
               x0(ix)=alog(alb)
               lx(ix)=1
             else
               print*,'Error in readapriori - alb must be > 0'
               stop
             endif
             sx(ix,ix) = (err/alb)**2
             print*,ix,err,alb,x0(ix),sx(ix,ix)
            enddo

            nx = nx+np


           elseif(varident(ivar,1).eq.889)then
C           **** Surface albedo scaling value *******
            jalb = nx+1
            ix = jalb
            read(27,*)alb,err
            if(alb.gt.0.0)then
               x0(ix)=alog(alb)
               lx(ix)=1
            else
               print*,'Error in readapriori - alb must be > 0'
               stop
            endif
            sx(ix,ix) = (err/alb)**2
            print*,ix,err,alb,x0(ix),sx(ix,ix)

            nx = nx+1

           elseif(varident(ivar,1).eq.777)then
C           **** Tangent altitude correction *******
            ix = nx+1
            read(27,*)x0(ix),err
            sx(ix,ix) = err**2
            jtan = ix

            nx = nx+1

           elseif(varident(ivar,1).eq.666)then
C           **** Pressure at given altitude
            ix = nx+1
            read(27,*)htan
            read(27,*)pre,err
            varparam(ivar,1) = htan
            if(pre.gt.0.0)then
              x0(ix)=alog(pre)
              lx(ix)=1
            else
              print*,'Error in readapriori - pre must be > 0'
              stop
            endif
            sx(ix,ix) = (err/pre)**2
            jpre = ix

            nx = nx+1

C **************** add radius variable : JM ***************
           elseif(varident(ivar,1).eq.555)then
C           **** Radius of planet *******
            ix = nx+1
            read(27,*)x0(ix),err
            sx(ix,ix) = err**2
            jrad = ix
          
            call readrefiplan(opfile,iplanet,xlat,radius)
            jradf = jrad
            radius2 = radius+x0(ix)

C	    print*,jrad,ix,sx(ix,ix),'jm2'

            nx = nx+1

           elseif (varident(ivar,1).eq.444)then
C            ** Variable cloud particle size distribution and composition
             read(27,1)rifile

             open(28,file=rifile,status='old')
C              Read mean radius and error
               read(28,*)r0,er0
               ix=nx+1
               x0(ix)=alog(r0)
               sx(ix,ix)=(er0/r0)**2
               lx(ix)=1
C              Read radius variance and error
               read(28,*)dr,edr
               ix=nx+2
               x0(ix)=alog(dr)
               sx(ix,ix)=(edr/dr)**2
               lx(ix)=1

C              Read number of pressures and correlation length (wavelength/
C				wavenumbers)
               read(28,*)np,clen
               varparam(ivar,1)=np


               call get_xsecA(opfile,nmode,nwave,wave,xsec)
               if(np.ne.nwave)then
       print*,'Error in readapriori.f. Number of wavelengths in ref.'
       print*,'index file does not match number of wavelengths in'
       print*,'xsc file. Wavelengths in these two files should match'
                print*,rifile
                print*,opfile
                stop
               endif

               varparam(ivar,2)=clen

C              read reference wavelength and nr at that wavelength
               read(28,*)vm,nm
               varparam(ivar,3)=vm
               varparam(ivar,4)=nm

C              read x-section normalising wavelength (-1 to not normalise)
               read(28,*)lambda0
               varparam(ivar,5)=lambda0
               do i=1,np             
                read(28,*)vi(i),nimag,err
                ix=nx+2+i
                x0(ix)=alog(nimag)
                sx(ix,ix)=(err/nimag)**2
                lx(ix)=1
               enddo
             close(28)

             do i=1,np
              do j=i+1,np
               delv = vi(i)-vi(j)
               arg = abs(delv/clen)
C               xfac = exp(-arg)
               xfac = exp(-arg*arg)
               if(xfac.ge.SXMINFAC)then  
                sx(nx+2+i,nx+2+j)=
     & sqrt(sx(nx+2+i,nx+2+i)*sx(nx+2+j,nx+2+j))*xfac
                sx(nx+2+j,nx+2+i)=sx(nx+2+i,nx+2+j)
               endif
              enddo
             enddo

             nx = nx+2+np
 
C **************** add mass variable  ***************
           elseif(varident(ivar,1).eq.333)then
C           **** surface ln(g) of planet *******
            ix = nx+1
            read(27,*)x0(ix),err
            sx(ix,ix) = err**2
            jlogg = ix
            jloggf = jlogg

            call readrefiplan(opfile,iplanet,xlat,radius)
            mass2 = 1e-20*10**(x0(ix))*(radius**2)/Grav

            nx = nx+1

           else
            print*,'vartype not recognised'
            stop
           endif

         endif

         if(nx.gt.mx)then 
           print*,'NX > MX ',nx,mx
           stop
         endif

10    continue
C     If both mass and radius being retrieved then we need to update
C     the mass using the a priori log(g) AND radius
      if(jrad.gt.0.and.jlogg.gt.0)then
         mass2 = 1e-20*10**(x0(jlogg))*(radius2**2)/Grav
      endif

      print*,'Total number of variables (nx) : ',nx

      close(27)


      if(lin.eq.2.or.lin.eq.3)then
       print*,'Readapriori: previous retrieval being used to
     1 update apriori'

       call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,
     1  varparamx,jsurfx,jalbx,jtanx,jprex,jradx,jloggx,nxx,xnx,
     2  sxx)

       xdiff = abs(xlat-xlatx)
       if(xdiff.gt.5.0)then
         print*,'Readapriori: Aborting - latitudes inconsistent'
         print*,xlatx,xlat
         stop
       endif
       if(nprox.ne.npro)then
         print*,'Readapriori: Aborting - npro inconsistent'
         print*,nprox,npro
         stop
       endif

       ioffx=0
       do 21 ivarx=1,nvarx
        npx=1
        if(varidentx(ivarx,1).le.100)then
          npx=npvar(varidentx(ivarx,3),npro)
        endif
        if(varidentx(ivarx,1).eq.888)npx=int(varparamx(ivarx,1))
        if(varidentx(ivarx,1).eq.444)npx=2+int(varparamx(ivarx,1))
     
        ioff=0
        do 22 ivar=1,nvar
         np=1
         if(varident(ivar,1).le.100)then
           np=npvar(varident(ivar,3),npro)
         endif
         if(varident(ivar,1).eq.888)np=int(varparam(ivar,1))
         if(varident(ivar,1).eq.444)np=2+int(varparam(ivar,1))


         if(varidentx(ivarx,1).eq.varident(ivar,1))then
          if(varidentx(ivarx,2).eq.varident(ivar,2))then
           if(varidentx(ivarx,3).eq.varident(ivar,3))then

            print*,'Updating variable : ',ivar,' :  ',
     1		(varident(ivar,j),j=1,3)
            do 33 i=1,np
             x0(ioff+i)=xnx(ioffx+i)
             do 34 j=1,np
              sx(ioff+i,ioff+j)=sxx(ioffx+i,ioffx+j)
34           continue
33          continue
            
           endif
          endif
         endif

         ioff=ioff+np

22      continue

        ioffx = ioffx+npx
 
21     continue

      endif

C     Write out x-data to temporary .str file for later routines.
      if(lin.eq.3)then
       call writextmp(runname,xlatx,nvarx,varidentx,varparamx,nprox,
     1  nxx,xnx,sxx,jsurfx,jalbx,jtanx,jprex,jradx,jloggx)
      endif

      return

      end
