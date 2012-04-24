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
C
C     *******************************************************************
      implicit none
      integer j,nmu
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer lowbc,liscat(maxcon),lncons(maxcon),lnorm(maxcon)
      integer isol,ncont,nf,i
      double precision mu(maxmu), wt(maxmu)
      real galb,dist,lcons(maxcon,maxscatpar),aphi,sol_ang,emiss_ang
      character*30 header
      character*30 scatfile(maxcon)
      character*100 radfile

      data scatfile/'hgphase1.dat','hgphase2.dat','hgphase3.dat',
     1	'hgphase4.dat','hgphase5.dat','hgphase6.dat','hgphase7.dat',
     2  'hgphase8.dat','hgphase9.dat','hgphase10.dat'/

1     format(a)
3     format(1X,a30)

      call file(radfile,radfile,'sca')
 

      do i=1,ncont
       liscat(i)=5
       lnorm(i)=1
       lncons(i)=0
       do j=1,10
        lcons(i,j)=0.0
       enddo
      end do


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
 



