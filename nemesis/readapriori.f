      subroutine readapriori(opfile,lin,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac,nx,x0,
     2  sx,lx)
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
C 	jxsc		integer		Position of start of x-section
C					spectrum
C	jtan		integer		Position of tangent altitude
C					correction
C	jpre		integer		Position of ref. tangent  pressure
C       jrad            integer         Position of radius of planet
C       jlogg           integer         Position of surface log_10(g) of planet
C       jfrac		integer		Position of fractional coverage
C	nx 		integer 	number of elements in state vector
C	x0(mx)		real		a priori vector
C	sx(mx,mx)	real		a priori covariance matrix
C	lx(mx)		integer		Log flag. 0 if real number
C						  1 if log number 
C	csx(mvar)	real		Ratio volume of shell/total volume of particle
C					for Maltmieser coated sphere model
C					For homogeneous sphere model, csx(ivar)=-1  
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
      integer jlogg,nmode,nwave,max_mode, max_wave,jxsc,icloud
      integer jfrac,jfracx
      parameter (max_mode = 10)
      parameter (max_wave = 1000)
      parameter(maxlat=100)

C     ****************************************************************
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/planrad.f'
      include '../radtran/includes/emcee.f'
      include 'arraylen.f'
C     ****************************************************************

      real xsec(max_mode,max_wave,2),wave(max_wave)
      real x0(mx),sx(mx,mx),err,err1,ref1,pref(maxpro)
      real eref(maxlat,maxpro),reflat(maxlat),htan
      real delp,xfac,pknee,eknee,edeep,xdeep,xlat,xlatx,xlonx,pshld
      real xstep,estep,efsh,xfsh,varparam(mvar,mparam),flat,hknee,pre
      real ref(maxlat,maxpro),clen,SXMINFAC,arg,valb,alb,refp,errp
      real xknee,xrh,erh,xcdeep,ecdeep,radius,Grav,plim
      real xcwid,ecwid,ptrop,refradius,xsc
      parameter (Grav=6.672E-11)
C     SXMINFAC is minimum off-diagonal factor allowed in the
C     a priori covariance matrix
      parameter (SXMINFAC = 0.001)
      real varparamx(mvar,mparam),xnx(mx),sxx(mx,mx),xdiff
      integer varident(mvar,3),ivar,nvar,nlevel,lin,jsurfx
      integer jalbx,jtanx,jprex,jradx,jlat,ilat,nlat,lx(mx)
      integer nprox,nvarx,varidentx(mvar,3),lpre,ioffx,ivarx
      integer npx,ioff,icond,npvar,jloggx,iplanet,jxscx,jlev
      character*100 opfile,buffer,ipfile,runname,rifile,xscfil
      integer nxx,nsec,ncont1,nlay,tmp
      real xwid,ewid,y,y0,lambda0,vi(mx),vtmp(mx),xt(mx)
      real r0,er0,dr,edr,vm,nm,nmshell,nimag,delv,xy
      real xldeep,eldeep,xlhigh,elhigh
      real v1,v1err,v2,v2err,p1,p1err,p2,p2err
      real tau0,ntemp,teff,alpha,T0,xf,exf
      real etau0,entemp,eteff,ealpha,eT0
      real csx,cserr,nimagshell,errshell,flon
      logical filexist
      integer i1,j1,nlocation,ilocation
      real tmpgrid(mparam),findgen(mparam),yout

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
      jxsc = -1
      jtan = -1
      jpre = -1
      jrad = -1
      jlogg = -1
      jfrac = -1
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
         csx = -1.0			!particles are assumed to be homogeneous by default
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

           elseif (varident(ivar,3).eq.-1)then
C          * continuous cloud, but cloud retrieved as particles/cm3 rather than 
C          * particles per gram to decouple it from pressure.
C          ********* continuous particles/cm3 profile ************************
             if(varident(ivar,1).ge.0)then
               print*,'readapriori.f error.'
               print*,'This model type is only for use with clouds'
               stop
             endif

             read(27,1)ipfile
             print*,'reading variable ',ivar,' from ',ipfile
             open(28,file=ipfile,status='old')
             read(28,*)nlevel,clen
             if(nlevel.ne.npro)then
              print*,'profiles must be listed on same grid as .prf'
              stop
             endif
             do 47 i=1,nlevel
              read(28,*)pref(i),ref(1,i),eref(1,i)
              ix = i+nx
C              **** cloud,take logs ***
               if(ref(1,i).gt.0.0) then
                 x0(ix) = alog(ref(1,i)) 
                 lx(ix)=1
               else 
                 print*,'Error in readapriori.f. Cant take log of zero'
                 print*,i,ref(i,i)
                 stop
               endif
               err = eref(1,i)/ref(1,i)
               sx(ix,ix)=err**2
47           continue
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
               delp = log(pref(j))-log(pref(i))
               
               arg = abs(delp/clen)
               xfac = exp(-arg)
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
             if(MCMCflag.eq.1)then
              if(varident(ivar,1).eq.-1)then
               hknee = MCMChknee
               xdeep = MCMCdeep
               xfsh = MCMCfsh
              else
               hknee = MCMChknee2
               xdeep = MCMCdeep2
               xfsh = MCMCfsh2
              endif

              eknee = 0.5*MCMChknee
              edeep = 0.5*MCMCdeep
              efsh = 0.5*MCMCfsh
             endif

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
C            ** profile held as peak OD, altitude and FWHM (km) (Gaussian) **
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

           elseif (varident(ivar,3).eq.17)then
C            ******** profile held as a value at a certain pressure and 
C            ******** fractional scale height, limited at pressures
C            ******** less than that specified 
C            Read in xknee,fsh,pknee,plim
             read(27,*)pknee,plim
             read(27,*)xknee,eknee
             read(27,*)xfsh,efsh
             varparam(ivar,1) = pknee
             varparam(ivar,2) = plim
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

           elseif (varident(ivar,3).eq.18)then
C            ******** profile held as a value at a certain pressure and 
C            ******** fractional scale height. vmr is constant at altitudes
C	     ******** above the reference level. 
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


           elseif (varident(ivar,3).eq.19)then
C            ******** cloud profile held as total optical depth plus
C            ******** base height and fractional scale height. Below the knee
C            ******** pressure the profile is set to zero - a simple
C            ******** cloud in other words! However, in this model, the 
C            ******** opacity is forced to zero at pressures  < 0.1 atm
C            ******** using a cut-off parameter
C            Read in xdeep,fsh,pknee
             read(27,*)hknee,eknee
             read(27,*)xdeep,edeep
             read(27,*)xfsh,efsh
             read(27,*)xcwid,ecwid

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

             ix = nx+4
             if(xcwid.gt.0.0)then
               x0(ix) = alog(xcwid)
               lx(ix)=1
             else
               print*,'Error in readapriori - xcwid must be > 0'
               stop
             endif
             sx(ix,ix) = (ecwid/xcwid)**2

             nx = nx+4

           elseif (varident(ivar,3).eq.20)then
C            ******** profile held as deep amount, fsh, knee pressure
C            ******** and tropopause cut-off pressure 
C            Read in xdeep,fsh,pknee,ptrop
             read(27,*)pknee,ptrop
             read(27,*)xdeep,edeep
             read(27,*)xfsh,efsh
             varparam(ivar,1) = pknee
             varparam(ivar,2) = ptrop
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

           elseif (varident(ivar,3).eq.21)then
C            ******** cloud profile held as total optical depth plus
C            ******** base height. Below the knee
C            ******** pressure the profile is set to zero
C            ******** Fractional scale height is scaled from that
C            ******** guessed by the radius of the associated size
C            ******** distribution in the 444/445 section.
C            Read in xdeep,fsh,pknee
             read(27,*)hknee,eknee
             read(27,*)xdeep,edeep
             read(27,*)xfsh,refradius
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
             x0(ix) = hknee
             sx(ix,ix) = eknee**2

             varparam(ivar,1)=xfsh
             varparam(ivar,2)=refradius

             nx = nx+2

           elseif (varident(ivar,3).eq.22)then
C           Parameterised Brown Dwarf T-profile
C            Read in tau0, n, Teff, alpha, T0
             read(27,*)tau0,etau0
             read(27,*)ntemp,entemp
             read(27,*)teff,eteff
             read(27,*)alpha,ealpha
             read(27,*)T0,eT0
             ix=nx+1
             x0(ix)=alog(tau0)
             lx(ix)=1
             err = etau0/tau0
             sx(ix,ix)=err**2
             ix=nx+2
             x0(ix)=alog(ntemp)
             lx(ix)=1
             err = entemp/ntemp
             sx(ix,ix)=err**2
             ix=nx+3
             x0(ix)=alog(teff)
             lx(ix)=1
             err = eteff/teff
             sx(ix,ix)=err**2
             ix=nx+4
             x0(ix)=alog(alpha)
             lx(ix)=1
             err = ealpha/alpha
             sx(ix,ix)=err**2
             ix=nx+5
             x0(ix)=alog(T0)
             lx(ix)=1
             err = eT0/T0
             sx(ix,ix)=err**2

             nx=nx+5

           elseif ( (varident(ivar,3).eq.23) .or. 
     >              (varident(ivar,3).eq.26) ) then
C            ******** 2 point gradient profile (NAT)
c		Profile is defined by two (p,v) points, with a linear gradient (in log p)
c		in between. The low pressure point is at (p1,v1) and the high pressure point 
c		is at (p2,v2). 
c		23: Profile is constant above/below this gradient region (i.e.
c		p<p1 v=v1 and p>p2 v=v2.) 
c		26: Profile is constant above this gradient region and zero below (i.e.
c		p<p1 v=v1 and p>p2 v=0.) 
c		All variable are retrieved. 
c		Not yet fully implemented for T	
             read(27,*) v1,v1err
             read(27,*) p1,p1err
             read(27,*) v2,v2err
             read(27,*) p2,p2err
             if (p1.gt.p2)then
               print*,'Error p1 should be less than p2'
               stop
             endif
c            ** param1 (low pressure limit)
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)= v1
              err   = v1err
             else
C             *** vmr, fcloud, para-H2 or cloud, take logs *********
              if(v1.gt.0.0)then
                x0(ix)=alog(v1)
                lx(ix)=1
              else
                print*,'Error in readapriori. v1 must be > 0.0'
                stop
              endif
              err = v1err/v1
             endif
             sx(ix,ix)=err**2
c            ** param2 (low pressure transition pressure)
             ix = nx+2
             if(p1.gt.0.0)then
                x0(ix)=alog(p1)
                lx(ix)=1
             else
                print*,'Error in readapriori. p1 must be > 0.0'
                stop
             endif
             err = p1err/p1
             sx(ix,ix)=err**2
c            ** param3 (high pressure limit)
             ix = nx+3
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)= v2
              err   = v2err
             else
C             *** vmr, fcloud, para-H2 or cloud, take logs *********
              if(v2.gt.0.0)then
                x0(ix)=alog(v2)
                lx(ix)=1
              else
                print*,'Error in readapriori. v2 must be > 0.0'
                stop
              endif
              err = v2err/v2
             endif
             sx(ix,ix)=err**2
c            ** param4 (high pressure transition pressure)
             ix = nx+4
             if(p2.gt.0.0)then
                x0(ix)=alog(p2)
                lx(ix)=1
             else
                print*,'Error in readapriori. p must be > 0.0'
                stop
             endif
             err = p2err/p2
             sx(ix,ix)=err**2
             
             nx = nx+4

           elseif (varident(ivar,3).eq.24)then
C            ******** profile held as variable knee pressure, one fixed value below knee and other fixed value above knee (REVIEW THIS IF DOING TEMP RETRIEVALS)
C            Read in xdeep,fsh,pknee
             read(27,*)pknee,eknee
             read(27,*)xdeep,edeep
             read(27,*)xstep,estep
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
             if(xstep.gt.0.0)then
               x0(ix) = alog(xstep)
               lx(ix)=1
             else
               print*,'Error in readapriori - xstep must be > 0'
               stop
             endif
             sx(ix,ix) = (estep/xstep)**2
             ix = nx+3
             x0(ix) = alog(pknee)
             lx(ix)=1

             sx(ix,ix) = (eknee/pknee)**2

             nx = nx+3

           elseif (varident(ivar,3).eq.25)then
C           Continuous profile but represented with fewer points than in .prf to achieve 
C           implicit smoothing and faster retrieval times
C            Read in number of points and any cross-correlation


             do i=1,np
C             References to pref can probably be removed entirely from this section,
C             but this is a temporary fix just to get rid of compiler issues
              pref(i)=0
             enddo
             read(27,1)ipfile
             print*,'reading variable ',ivar,' from ',ipfile
             open(28,file=ipfile,status='old')
C            np = number of points over which to retrieve profile (as opposed to npro which is total number of points in .prf)
C            nlay = number of homogeneous layers in .drv file
             read(28,*)np,clen,nlay
             varparam(ivar,1)=np
             varparam(ivar,2)=nlay
             if(np.gt.npro)then
              print*,'Error readapriori:'
              print*,'np > npro'
              stop
             endif
             do i=1,np
              read(28,*)pknee,xdeep,edeep!pressure, aerosol density, error
              ix = i+nx
C             For vmrs and cloud density always hold the log. 
C             Avoids instabilities arising from greatly different profiles 
C             such as T and vmrs
              if(varident(ivar,1).eq.0)then
C              *** temperature, leave alone ****
               x0(ix)=xdeep
               sx(ix,ix)=edeep**2
              else
C              **** vmr, cloud, para-H2 , fcloud, take logs ***
               if(pknee.lt.0.0) then
                 print*,'Error in readapriori.f.'
                 print*,'Cannot have negative pressures'
                 stop
               endif
               if(xdeep.lt.0.0) then
                 print*,'Error in readapriori.f.'
                 print*,'Cannot take log of negative values'
                 stop
               endif
               x0(ix)=alog(xdeep)
               sx(ix,ix)=(edeep/xdeep)**2
               lx(ix)=1
              endif
              varparam(ivar,i+2)=pknee

             enddo
             close(28)


c             print*,'pref = ', pref
c             stop

             do i = 1,np
              ix = nx + i         
              do j = i+1,np
               jx = nx + j
               if(pref(i).lt.0.0) then
                 print*,'Error in readapriori.f. A priori file '
                 print*,'must be on pressure grid '
                 stop
               endif
               if(np.gt.mparam-2)then
                print*,'Error readapriori.f: np + 2 > MPARAM'
                print*,'Either reduce value of np in ',ipfile
                print*,'or increase value of MPARAM (arraylen.f)'
                stop
               endif            
c               print*,'DIAG1',pref(i),pref(j),alog(pref(j))
               delp = log(varparam(ivar,j+2))-log(varparam(ivar,i+2))
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

             nx=nx+np

C            Write out pressure.lay for layer splitting in .drv file according to pressure grid specified in this parametrisation
C            Note that if more than one VARIDENT(IVAR,3) = 25 profile to be retrieved, pressure.lay will only be computed for the first profile listed in jupiter.apr
C            Therefore if more than one VARIDENT(IVAR,3) = 25 profile to be retrieved, best to make the first profile the one with the highest pressure resolution
             inquire(file='pressure.lay',exist=filexist)
             if(filexist.eqv..FALSE.)then
              open(901,file='pressure.lay',status='unknown')
               write(901,*)'*********Header***********'
               write(901,*)nlay
               if(mod(nlay,np-1).eq.0)then
                do i=3,np+1
                 delp = log(varparam(ivar,i+1))-log(varparam(ivar,i))
                 write(901,*)varparam(ivar,i)
                 tmp = nlay/(np-1)
                 do j=1,tmp-1
                  write(901,*)exp(log(varparam(ivar,i)) + j*delp/tmp)
                 enddo
                enddo
               else
C               Interpolate homogeneous levels from pressure grid
                do i=1,np
                 tmpgrid(i) = log(varparam(ivar,i+2))
                 findgen(i) = float(i-1)
                enddo
                do i=1,nlay
                 call verint(findgen,tmpgrid,np,yout,
     1                        (i-1)*(findgen(np)/float(nlay)))
                 print*,i,(i-1)*(findgen(np)/float(nlay)),exp(yout)
                 write(901,*)exp(yout)
                enddo
               endif
              close(901)
             endif

           elseif (varident(ivar,3).eq.27) then
C            step profile
             read(27,*)pknee,eknee
             read(27,*)xdeep,edeep
             read(27,*)xstep,estep

c		 xn1=log(deep)		[=unlogged deep value if temperature]
c		 xn2=log(shallow)		[=unlogged shallow value if temperature]
c		 xn3=log(knee pressure)             

c		 ** Deep value **
             ix = nx+1
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix)= xdeep
              err   = edeep
             else
C             *** vmr take logs *********
              if(xdeep.gt.0.0)then
                x0(ix) = alog(xdeep)
                lx(ix) = 1
              else
                print*,'Error in readapriori. xdeep must be > 0.0'
                stop
              endif
              err = edeep/xdeep
             endif
             sx(ix,ix)=err**2
             
c		 ** Shallow value **
             ix = nx+2
             if(varident(ivar,1).eq.0)then
C             *** temperature, leave alone ********
              x0(ix) = xstep
              err    = estep
             else
C             *** vmr take logs *********
              if(xstep.gt.0.0)then
                x0(ix) = alog(xstep)
                lx(ix) = 1
              else
                print*,'Error in readapriori. xstep must be > 0.0'
                stop
              endif
              err = estep/xstep
             endif
             sx(ix,ix)=err**2
             
C		 ** Knee **
             ix = nx+3
             x0(ix) = alog(pknee)
             lx(ix) = 1
             sx(ix,ix) = (eknee/pknee)**2

             nx = nx+3

           elseif(varident(ivar,3).eq.28) then
C          Modify one element of a profile

            ix = nx+1
            read(27,*)jlev
            varparam(ivar,1) = jlev
            read(27,*)refp,errp

            if(varident(ivar,1).eq.0)then
c             *** temperature, leave alone ***
              x0(ix) = refp
              err = errp 
            else
c             *** vmr, cloud, para-H2 , fcloud, take logs ***
              if(refp.gt.0.0)then
                x0(ix) = alog(refp)
                lx(ix)=1
              else
                 print*,'Error in readapriori.f. Cant take log of zero'
                 print*,refp
                 stop
              endif
              err = errp/refp
            endif
    
            sx(ix,ix) = err**2
            nx = nx+1


           elseif(varident(ivar,3).eq.29)then
C          ***** multiple continuous profiles at different locations  ********
             read(27,1)ipfile
             print*,'reading variable ',ivar,' from ',ipfile
             open(28,file=ipfile,status='old')
             read(28,*)nlocation,nlevel,clen
             if(nlevel.ne.npro)then
              print*,'profiles must be listed on same grid as .prf'
              stop
             endif
             if(2*nlocation+1.gt.mparam)then
              print*,'2*nlocation+1 > mparam',nlocation,mparam
              print*,'Need to reduce nlocation or increase mparam'
              stop
             endif
             varparam(ivar,1)=nlocation
             i=2
             do 201 ilocation=1,nlocation
C             Read in latitude/longitude
              read(28,*)flat,flon
              varparam(ivar,i)=flat
              varparam(ivar,i+1)=flon
              i=i+2
201          continue

             do 202 ilocation=1,nlocation
C             Read in apriori fraction
              do 203 i=1,nlevel
               read(28,*)pref(i),ref(1,i),eref(1,i)
               ix = i+nx
C              For vmrs and cloud density always hold the log. 
C              Avoids instabilities arising from greatly different profiles 
C              such as T and vmrs
               if(varident(ivar,1).eq.0)then
C               *** temperature, leave alone ****
                x0(ix) = ref(1,i)
                err = eref(1,i)
               else
C               **** vmr, cloud, para-H2 , fcloud, take logs ***
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

203           continue

              do i = 1,nlevel
               ix = nx + i 
               do j = i+1,nlevel
                jx = nx + j
                if(pref(i).lt.0.0) then
                  print*,'Error in readapriori.f. A priori file '
                  print*,'must be on pressure grid '
                  stop
                endif            

                delp = log(pref(j))-log(pref(i))
                
                arg = abs(delp/clen)
                xfac = exp(-arg)
C                xfac = exp(-arg*arg)
                if(xfac.ge.SXMINFAC)then  
                 sx(ix,jx) = sqrt(sx(ix,ix)*sx(jx,jx))*xfac
                 sx(jx,ix) = sx(ix,jx)
                endif
               enddo
              enddo

              nx=nx+nlevel

202          continue

C             nx=nx+nlevel

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
            read(27,*)np,clen
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

            do i=1,np
             do j=1,np
               delv = vi(i)-vi(j)
               arg = abs(delv/clen)
               xfac = exp(-arg)
               if(xfac.ge.SXMINFAC)then  
                sx(nx+2+i,nx+2+j)=
     & sqrt(sx(nx+2+i,nx+2+i)*sx(nx+2+j,nx+2+j))*xfac
                sx(nx+2+j,nx+2+i)=sx(nx+2+i,nx+2+j)
               endif
             enddo
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

           elseif(varident(ivar,1).eq.887)then
C           **** Cloud x-section spectrum *******
C           Read in number of points, cloud id, and correlation between elements.
            read(27,*)np, icloud, clen
            varparam(ivar,1)=np
            varparam(ivar,2)=icloud
            jxsc = nx+1
            do i=1,np             
             ix = jxsc+i-1
             read(27,*)vi(i),xsc,err
             if(xsc.gt.0.0)then
               x0(ix)=alog(xsc)
               lx(ix)=1
             else
               print*,'Error in readapriori - xsc must be > 0'
               stop
             endif
             sx(ix,ix) = (err/xsc)**2
             print*,ix,err,xsc,x0(ix),sx(ix,ix)
            enddo

C           Check the wavelengths are consistent with the .rxs file
            call file(runname,xscfil,'rxs')
            open(28,file=xscfil,status='old')
54          read(28,1)buffer
            if(buffer(1:1).eq.'#')goto 54
            read(buffer,*)ncont1
            j=0
105         j=j+1
            read(28,*,end=106)vtmp(j),(xt(i),i=1,ncont1)
            read(28,*,end=106)(xt(i),i=1,ncont1)
            xy=abs(vi(j)-vtmp(j))
            if(xy.gt.0.01) then 
             print*,'Possible problem in readapriori.f. Model 887'
             print*,'Wavelngth/wavenumb. in .rxs inconsistent with .apr'
             print*,j,vi(j),vtmp(j)
            endif
            goto 105
106         continue
            nsec=j-1
            if (nsec.ne.np)then
             print*,'Error in readapriori.f'
             print*,'Model 887: np = ',np
             print*,'nsec in .rxs file = ',nsec
             stop
            endif
            close(28)

            do i=1,np
             do j=1,np
               delv = vi(i)-vi(j)
               arg = abs(delv/clen)
               xfac = exp(-arg)
               if(xfac.ge.SXMINFAC)then  
                sx(nx+2+i,nx+2+j)=
     & sqrt(sx(nx+2+i,nx+2+i)*sx(nx+2+j,nx+2+j))*xfac
                sx(nx+2+j,nx+2+i)=sx(nx+2+i,nx+2+j)
               endif
             enddo
            enddo
   
            nx = nx+np

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
            print*,'radius2 = ',radius2
            nx = nx+1

           elseif (varident(ivar,1).eq.444)then
C            ** Variable cloud particle size distribution and composition
             read(27,1)rifile

             open(28,file=rifile,status='old')
C              Read mean radius and error
               read(28,*)r0,er0
               ix=nx+1

               if(MCMCflag.eq.1)then
                r0 = MCMCpr
                er0 = 1.0e-7*r0
               endif

               x0(ix)=alog(r0)
               sx(ix,ix)=(er0/r0)**2
               lx(ix)=1
C              Read radius variance and error
               read(28,*)dr,edr

               if(MCMCflag.eq.1)then
                dr = MCMCpvar
                edr = 1.0e-7*dr
               endif

               ix=nx+2
               x0(ix)=alog(dr)
               sx(ix,ix)=(edr/dr)**2
               lx(ix)=1

C              Read number of wavelengths and correlation length (wavelength/
C				wavenumbers)
               read(28,*)np,clen
               varparam(ivar,1)=np


               call get_xsecA(opfile,nmode,nwave,wave,xsec)
               if(np.ne.nwave)then
       print*,'Error in readapriori.f. Number of wavelengths in ref.'
       print*,'index file does not match number of wavelengths in'
       print*,'xsc file. Wavelengths in these two files should match'
                print*,rifile,np
                print*,opfile,nwave
                stop
               endif

               varparam(ivar,2)=clen

C              read reference wavelength and nr at that wavelength
               read(28,*)vm,nm
               if(MCMCflag.eq.1)then
                nm = MCMCreal
               endif
               varparam(ivar,3)=vm
               varparam(ivar,4)=nm

C              read x-section normalising wavelength (-1 to not normalise)
               read(28,*)lambda0
               varparam(ivar,5)=lambda0
               do i=1,np             
                read(28,*)vi(i),nimag,err
                if(MCMCflag.eq.1)then
                 nimag = MCMCimag(i)
                 err = 0.1*nimag
                endif
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
               xfac = exp(-arg)
C               xfac = exp(-arg*arg)
               if(xfac.ge.SXMINFAC)then  
                sx(nx+2+i,nx+2+j)=
     & sqrt(sx(nx+2+i,nx+2+i)*sx(nx+2+j,nx+2+j))*xfac
                sx(nx+2+j,nx+2+i)=sx(nx+2+i,nx+2+j)
               endif
              enddo
             enddo

             nx = nx+2+np

           elseif (varident(ivar,1).eq.445)then
C            ** Variable cloud particle size distribution and composition using Maltmieser coated sphere model **
             read(27,1)rifile

             open(28,file=rifile,status='old')	!open cloud.dat
C              Read mean radius and error
               read(28,*)r0,er0
               ix=nx+1				!nx = number of elements in state vector
               x0(ix)=alog(r0)			!set a priori vector at position ix to the log of the mean radius
               sx(ix,ix)=(er0/r0)**2		!set a priori covar matrix at position (ix,ix) to (error/mean radius)^2
               lx(ix)=1				!set log flag
C              Read radius variance and error
               read(28,*)dr,edr
               ix=nx+2
               x0(ix)=alog(dr)			!do the same with the radius variance at the next position in the a priori vector and covar matrix
               sx(ix,ix)=(edr/dr)**2
               lx(ix)=1

C              Read ratio of shell volume wrt total volume of particle, with error                
               read(28,*)csx,cserr
               if((csx.ge.1).or.(csx.le.0))then
                print*,'Error readapriori:'
                print*,'csx must be between 0 and 1'
                print*,'csx = ', csx
                stop
               endif
               if((csx+cserr.ge.1).and.(csx-cserr.le.0))then
                print*,'Error readapriori:'
                print*,'cserr too large'
                print*,'csx+cserr = ', csx+cserr
                print*,'csx-cserr = ', csx-cserr
                stop
               endif
               ix=nx+3
               x0(ix)=alog(csx)
               sx(ix,ix)=(cserr/csx)**2
               lx(ix)=1

C              Read number of wavelengths and correlation length (wavelength/
C				wavenumbers)
               read(28,*)np,clen
               varparam(ivar,1)=np

C		Read .xsc file and get extinction cross-sec and sing-scat albedo for each aerosol type (both output in variable xsec) for each wavelength (wave)
               call get_xsecA(opfile,nmode,nwave,wave,xsec)
               if(np.ne.nwave)then
       print*,'Error in readapriori.f. Number of wavelengths in ref.'
       print*,'index file does not match number of wavelengths in'
       print*,'xsc file. Wavelengths in these two files should match'
                print*,rifile,np
                print*,opfile,nwave
                stop
               endif

               varparam(ivar,2)=clen

C              read reference wavelength and nr at that wavelength (core,then shell)
               read(28,*)vm,nm,nmshell
               varparam(ivar,3)=vm
               varparam(ivar,4)=nm
               varparam(ivar,6)=nmshell

C              read x-section normalising wavelength (-1 to not normalise)
               read(28,*)lambda0
               varparam(ivar,5)=lambda0
               do i=1,np             
                read(28,*)vi(i),nimag,err,nimagshell,errshell	!wavelength, core imag RI, shell imag RI
                ix=nx+3+i
                x0(ix)=alog(nimag)		!set next position in a priori vector to core imag RI
                sx(ix,ix)=(err/nimag)**2	!set a priori covariance matrix
                lx(ix)=1
                x0(ix+np)=alog(nimagshell)  !ditto with shell imag RI
                sx(ix+np,ix+np)=(errshell/nimagshell)**2
                lx(ix+np)=1
               enddo
             close(28)

             do i=1,np
              do j=i+1,np
               delv = vi(i)-vi(j)
               arg = abs(delv/clen)
               xfac = exp(-arg)
C               xfac = exp(-arg*arg)
               if(xfac.ge.SXMINFAC)then  
                sx(nx+3+i,nx+3+j)=
     & sqrt(sx(nx+3+i,nx+3+i)*sx(nx+3+j,nx+3+j))*xfac
                sx(nx+3+j,nx+3+i)=sx(nx+3+i,nx+3+j)
                sx(nx+3+i+np,nx+3+j+np)=
     & sqrt(sx(nx+3+i+np,nx+3+i+np)*sx(nx+3+j+np,nx+3+j+np))*xfac
                sx(nx+3+j+np,nx+3+i+np)=sx(nx+3+i+np,nx+3+j+np)
               endif
              enddo
             enddo

             nx = nx+3+(np*2)

 
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
            print*,'mass2 = ',mass2
            nx = nx+1



C **************** Sromovsky Cloud Model  ***************
           elseif(varident(ivar,1).eq.222)then
C           **** Larry's discrete cloud model *******
C           LTC base pressure
            ix = nx+1
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           LTC opacity
            ix = nx+2
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC base pressure
            ix = nx+3
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC opacity
            ix = nx+4
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC base pressure
            ix = nx+5
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC opacity
            ix = nx+6
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Tropospheric haze opacity
            ix = nx+7
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Stratospheric haze opacity
            ix = nx+8
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           read in haze p1,p2,p3,n1,n2
            read(27,*)(varparam(ivar,j),j=1,5)

            nx = nx+8

C **************** Sromovsky Cloud Model with methane  ***************
           elseif(varident(ivar,1).eq.223)then
C           **** Larry's discrete cloud model with methane *******
C           LTC base pressure
            ix = nx+1
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           LTC opacity
            ix = nx+2
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC base pressure
            ix = nx+3
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC opacity
            ix = nx+4
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Methane depletion base pressure
            ix = nx+5
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Methane depletion factor
            ix = nx+6
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC opacity
            ix = nx+7
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Tropospheric haze opacity
            ix = nx+8
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Stratospheric haze opacity
            ix = nx+9
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           read in haze p1,p2,p3,n1,n2
            read(27,*)(varparam(ivar,j),j=1,5)

            nx = nx+8

C **************** Sromovsky Cloud Model  ***************
           elseif(varident(ivar,1).eq.224)then
C           **** Larry's discrete cloud model, but extendable UTC  *******
C           LTC base pressure
            ix = nx+1
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           LTC opacity
            ix = nx+2
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC base pressure
            ix = nx+3
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC opacity
            ix = nx+4
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC base pressure
            ix = nx+5
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC opacity
            ix = nx+6
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC fsh
            ix = nx+7
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Tropospheric haze opacity
            ix = nx+8
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Stratospheric haze opacity
            ix = nx+9
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           read in haze p1,p2,p3,n1,n2
            read(27,*)(varparam(ivar,j),j=1,5)

            nx = nx+9

C **************** Sromovsky Cloud Model  ***************
           elseif(varident(ivar,1).eq.225)then
C           **** Larry's discrete cloud model, but extendable UTC and
C           depleted methane  *******
C           LTC base pressure
            ix = nx+1
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           LTC opacity
            ix = nx+2
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC base pressure
            ix = nx+3
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           MTC opacity
            ix = nx+4
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC base pressure
            ix = nx+5
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC opacity
            ix = nx+6
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC fsh
            ix = nx+7
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC tropopause cut-off rate
            ix = nx+8
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Tropospheric haze opacity
            ix = nx+9
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Stratospheric haze opacity
            ix = nx+10
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           Methane depletion factor
            ix = nx+11
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           read in haze p1,p2,p3,n1,n2
            read(27,*)(varparam(ivar,j),j=1,5)

            nx = nx+11

C **************** Two Cloud Model  ***************
           elseif(varident(ivar,1).eq.226)then
C           **** Simple two cloud model *******

C           LTC base pressure
            ix = nx+1
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           LTC top pressure
            ix = nx+2
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           LTC fsh
            ix = nx+3
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           LTC opacity
            ix = nx+4
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1


C           UTC base pressure
            ix = nx+5
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC top pressure
            ix = nx+6
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC fsh
            ix = nx+7
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           UTC opacity
            ix = nx+8
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

            nx = nx+8

C********************** Creme Brulee Model *******************************
           elseif(varident(ivar,1).eq.227)then
C           TC base pressure
            ix = nx+1
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           TC opacity
            ix = nx+2
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           TC fsh
            ix = nx+3
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           TC top pressure/CB base pressure
            ix = nx+4
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           CB opacity
            ix = nx+5
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           SH base pressure
            ix = nx+6
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           SH opacity
            ix = nx+7
            read(27,*)r0,err
            x0(ix)=alog(r0)
            sx(ix,ix) = (err/r0)**2
            lx(ix)=1

C           CB top pressure, specified as a ratio of the CB base pressure.
C           (for instance, if set to 0.2 then CB top pressure = 0.2 x CB base pressure)
C           If set to 0, default Creme Brulee model is used where CB top pressure = 0.9 x CB base pressure)
	    read(27,*)varparam(ivar,1)
            if(varparam(ivar,1).gt.0.9)then
             print*,'Error readapriori:'
             print*,'Varparam(ivar,1) must be less than 0.9'
             stop
            endif

            nx = nx+7
           elseif(varident(ivar,1).eq.102)then
            ix = nx+1
            read(27,*)x0(ix),err
            sx(ix,ix) = err**2
            jfrac = ix
            nx=nx+1
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


      if(lin.eq.2.or.lin.eq.3.or.lin.eq.4)then
       print*,'Readapriori: previous retrieval being used to
     1 update apriori'

       call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,
     1  varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx,
     2  jfracx,nxx,xnx,sxx)

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
          npx=npvar(varidentx(ivarx,3),npro,varparamx(ivarx,1))
        endif
        if(varidentx(ivarx,1).eq.888)npx=int(varparamx(ivarx,1))
        if(varidentx(ivarx,1).eq.887)npx=int(varparamx(ivarx,1))
        if(varidentx(ivarx,1).eq.444)npx=2+int(varparamx(ivarx,1))
        if(varidentx(ivarx,1).eq.445)npx=3+int(varparamx(ivarx,1))
     
        ioff=0
        do 22 ivar=1,nvar
         np=1
         if(varident(ivar,1).le.100)then
           np=npvar(varident(ivar,3),npro,varparam(ivar,1))
         endif
         if(varident(ivar,1).eq.888)np=int(varparam(ivar,1))
         if(varident(ivar,1).eq.887)np=int(varparam(ivar,1))
         if(varident(ivar,1).eq.444)np=2+int(varparam(ivar,1))
         if(varident(ivar,1).eq.445)np=3+int(varparam(ivar,1))


         if(varidentx(ivarx,1).eq.varident(ivar,1))then
          if(varidentx(ivarx,2).eq.varident(ivar,2))then
           if(varidentx(ivarx,3).eq.varident(ivar,3))then

            if(varidentx(ivarx,3).eq.28)then
             if(varparamx(ivarx,1).eq.varparam(ivar,1))then
 
              print*,'Updating variable : ',ivar,' :  ',
     1           (varident(ivar,j),j=1,3)
              do i=1,np
               x0(ioff+i)=xnx(ioffx+i)
               do j=1,np
                sx(ioff+i,ioff+j)=sxx(ioffx+i,ioffx+j)
               enddo
              enddo

             endif             
            else

             print*,'Updating variable : ',ivar,' :  ',
     1 		(varident(ivar,j),j=1,3)
             do 33 i=1,np
              x0(ioff+i)=xnx(ioffx+i)
              do 34 j=1,np
               sx(ioff+i,ioff+j)=sxx(ioffx+i,ioffx+j)
34            continue
33           continue
            
            endif

           endif
          endif
         endif

         ioff=ioff+np

22      continue

        ioffx = ioffx+npx
 
21     continue

      endif

C     Write out x-data to temporary .str file for later routines.
      if(lin.eq.3.or.lin.eq.4)then
       call writextmp(runname,xlatx,nvarx,varidentx,varparamx,nprox,
     1  nxx,xnx,sxx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx,
     2  jfracx)
      endif

      return

      end
