      subroutine mod_scatter(radfile,ncont,nmu,mu,wt,isol,
     1dist,lowbc,galb,nf,sol_ang,emiss_ang,aphi)
C     $Id:
C     *******************************************************************
C     Subroutine to write the .sca file for a scattering CIRSradg run.
C
C     Input variables
C	radfile		character*100	run name
C	nmu		integer		Number of zenith ordinates
C	ncont		integer		Number of aerosol types
C	mu(maxmu)	double prec	Cos(Zenith Angles)
C	wt(maxmu)	double prec	Quadrature weights for zenith angles
C	isol		integer		Sun on/off
C	dist		real		Distance to Sun (AU)
C	lowbc		integer		Lower boundary condition
C	galb		real		Ground albedo
C	nf		integer		Number of Fourier azimuth coeffs.
C	sol_ang		real		Solar zenith angle
C	emiss_ang	real	Viewing zenith angle
C	aphi  		real	Azimuth angle
C    
C     Pat Irwin 4/4/01		Original
C     Pat Irwin 17/10/03	Tidied for Nemesis
c     20nov24     NAT: added options that can be commented out
c                 for MAXCON/MAXSCATPAR=4,14 (kept the 10 version too)
c                 Not satisfying to have this hard wired but I doubt it
c                 Needs to be changed regularly.
c                 NB MAXCON & MAXSCATPAR are set in includes/arrdef.f
c                 and must be the same. If they are changed there are
c                 Hardwired bits in scatter/phase1.f and nemesis/mod_scatter.f
c                 that also need changing. Currently 10,14 have options that 
c                 can be commented/uncommented as needed. 10 is default/original
C
C     *******************************************************************
      implicit none
      integer j,nmu
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer lowbc,liscat(maxcon),lncons(maxcon),lnorm(maxcon)
      integer isol,ncont,nf,i,imie
      double precision mu(maxmu), wt(maxmu)
      real galb,dist,lcons(maxcon,maxscatpar),aphi,sol_ang,emiss_ang
      character*30 header
      character*30 scatfile(maxcon),scatfile1(maxcon),scatfile0(maxcon)
      character*100 radfile
      common /imiescat/imie

c ** version for MAXCON=10
      data scatfile0/'hgphase1.dat','hgphase2.dat','hgphase3.dat',
     &	'hgphase4.dat','hgphase5.dat','hgphase6.dat','hgphase7.dat',
     &  'hgphase8.dat','hgphase9.dat','hgphase10.dat'/
      data scatfile1/'PHASE1.DAT','PHASE2.DAT','PHASE3.DAT',
     & 'PHASE4.DAT','PHASE5.DAT','PHASE6.DAT','PHASE7.DAT',
     & 'PHASE8.DAT','PHASE9.DAT','PHASE10.DAT'/
c ** version for MAXCON=14
c      data scatfile0/'hgphase1.dat','hgphase2.dat','hgphase3.dat',
c     &	'hgphase4.dat','hgphase5.dat','hgphase6.dat','hgphase7.dat',
c     &  'hgphase8.dat','hgphase9.dat','hgphase10.dat','hgphase11.dat',
c     &  'hgphase12.dat','hgphase13.dat','hgphase14.dat'/
c      data scatfile1/'PHASE1.DAT','PHASE2.DAT','PHASE3.DAT',
c     & 'PHASE4.DAT','PHASE5.DAT','PHASE6.DAT','PHASE7.DAT',
c     & 'PHASE8.DAT','PHASE9.DAT','PHASE10.DAT','PHASE11.DAT',
c     & 'PHASE12.DAT','PHASE13.DAT','PHASE14.DAT'/

1     format(a)
3     format(1X,a30)

      call file(radfile,radfile,'sca')
 
      if(imie.eq.1)then

       do i=1,ncont
        liscat(i)=4
        lnorm(i)=1
        lncons(i)=0
        scatfile(i)=scatfile1(i)
        do j=1,10
         lcons(i,j)=0.0
        enddo
       enddo

      else

       do i=1,ncont
        liscat(i)=5
        lnorm(i)=1
        lncons(i)=0
        scatfile(i)=scatfile0(i)
        do j=1,10
        lcons(i,j)=0.0
        enddo
       enddo

      endif

      open(49,file=radfile,status='unknown')
      header = 'Output from NEMESIS mod_scatter'
      write(49,1)header
      write(49,*)nmu,'                 ! Number of zenith ordinates'
      do 15 i=1,nmu
       write(49,*)mu(i),wt(i)
15    continue
      write(49,*)isol,'                           ! Sunlight switch'
      write(49,*)dist,'                       ! Solar distance (AU)'
      write(49,*)lowbc,'                 ! Lower boundary condition'
      write(49,*)galb,'                             ! Ground albedo'
      write(49,*)sol_ang,emiss_ang,'    ! Solar/emiss zenith angles'
      write(49,*)aphi,'                             ! Azimuth angle'
      write(49,*)nf,'                ! Number of fourier components'
      write(49,*)ncont,'                 ! Number of particle types'
      do 150 i=1,ncont
       write(49,*)liscat(i),'                             ! Scat ID'
       write(49,*)lnorm(i),'                              ! Norm ID'
       if(liscat(i).lt.4)then
        write(49,*)lncons(i),'                              ! ncons'
        do 35 j=1,lncons(i)
         write(49,*)lcons(i,j),'                             ! cons'
35      continue
       else
        write(49,3)scatfile(i)
       end if
150   continue

      close(49)
   

      return
      end
 



