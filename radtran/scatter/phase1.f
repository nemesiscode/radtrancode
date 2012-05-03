      subroutine phase1( calpha, p, iscat, cons, icons, icont, ncont, 
     1   vwave)
C     $Id: phase1.f,v 1.3 2011-06-17 15:57:54 irwin Exp $
C **********************************************************************
C
C     Subroutine to calculate the phase function P at a single phase angle 
C     CALPHA.
C
C     P is normalized to 1 over all solid angles.
C
C     Input variables:
C	calpha	real	cosine of scattering angle
C  	iscat	integer	0:  dipole		  -- 0 coefficients in CONS
C          		1:  isotropic scattering  -- 0 coefficients in CONS
C          		2:  Henyey-Greenstein     -- 3 coefficients in CONS
C	   		3:  Legendre polynomial   -- ICONS (<= 10) coefficients 
C						     in CONS
C	   		4:  read in from disk file 
C			5:  read in  Henyey-Greenstein from file
C	cons(maxscatpar)	real*8	Required constants for phase calc.
C	icons	integer	Number of elements of cons used.
C	icont	integer	Aerosol mode ID (for external files)
C	ncont	integer	Number of aerosol types
C	vwave	real	Calculation wavenumber (for external files)
C   
C     Output variable:
C	p	real	Phase function
C
C **********************************************************************
C
C_HIST:	7dec88 ... initial version, for H-G. Andy Collard and L.Kamp
C   	9dec88 ... added ISCAT with isotropic/dipole options
C  	21dec88 ... revised to use VENUS convention for ISCAT
C   	5jan89 ... fixed code for ISCAT
C   	9jan89 ... use value of PI identical to CLOUD's
C   	7feb89 ... added Legendre option (replaced subr. PHASLEG)
C   	8feb89 ... added option for disk file input (ISCAT=4),
C		   consolidated data parameters into CONS.
C   	1mar89 ... added ISCAT=5 option;  change to IMPLICIT NONE for
C		   use by pgm. VENUS;  open disk file in calling routine
C   	2mar89 ... ensure that first Legendre constant is 1/(4*pi) to full
C		   double precision accuracy
C  	25mar89 ... fixed bug re-initializing files for ISCAT=5
C	14oct02	PDP	Removed references to "exit" as the command is not
C			supported under the Intel FORTRAN Compilter
C			framework. Left the old code commented out just in
C			case "CALL EXIT" is not equivalent to "STOP".
C
C **********************************************************************

	IMPLICIT NONE
        include '../includes/arrdef.f'
        integer maxpts, ncont, icont, oldiwave, oldicont,idum1
        integer irec,lunit,j,iwave
	parameter (maxpts=100)
	real*4 xmu(maxpts), pfunc(maxpts),calph1 
        real*8 p, calpha, f1, f2, pi,
     1   hg11, hg12, hg21, hg22, xf, x0, x1, pa0, pa1, theta(maxpts)
	real*4 frq, oldfrq, dum1,dum2,dum3,vwave,dummy,p1,oldvwave
	integer icons, n0, n1, imu, npts, mod, oldmod, iscat, i, k
	logical init/.true./, initpi/.true./
	character*1 head
        character*512 buffer
	save pi, x0, x1, pa0, pa1, imu, initpi, init, 
     1  oldvwave,pfunc,xmu

        real*4 thet(50),phas(10,3,50),header(10,3,3),cthet(50)
        real*4 v1,v0,dv,vv,frac
        integer irec1(10),npoint,nphas,ncons,maxrec(10)

        common /phase/thet,cthet,phas,header,irec1,nphas,maxrec

C This dimension assumes that Legendre option needs the most constants:
	real*8 cons(maxscatpar)

C Coefficients & scale factors of Legendre polynomials P(n-1): (1st index
C of JCOEF is power of MU, second is N=0,1,2,...,MAXSCATPAR. Note that a
C polynomial of even N has only even powers, and odd N has only odd
C powers.)
	integer jcoef(maxscatpar/2,maxscatpar)/
	1 1, 4*0,
	2 1, 4*0,
	3 -1, 3, 3*0,
	4 -3, 5, 3*0,
	5 3, -30, 35, 2*0,
	6 15, -70, 63, 2*0,
	7 -5, 105, -315, 231, 0, 
	8 -35, 315, -693, 429, 0, 
	9 35, -1260, 6930, -12012, 6435, 
	1 315, -4620, 18018, -25740, 12155/
C Each scale factors is for polynomial 2i-1 & 2i:
	integer jdiv(maxscatpar/2)/ 1, 2, 8, 16, 128/


C **********************************************************************
	if (initpi) then
	  pi = 4.0d0*datan(1.0d0)
	  initpi = .false.
	endif

	calpha = dmin1( dmax1(calpha,-1.0d0), 1.0d0)

	goto (10,20,30,40,50,60),iscat+1
	PRINT*,' PHASE1.f :: Error invalid scattering option.'
	PRINT*,' PHASE1.f :: Stopping program.'
	STOP
cc	print*,' PHASE1:  invalid scattering option'
cc	call exit

C **********************************************************************
C ISCAT= 0 : Dipole
10	p = 0.75d0*(1.0d0+calpha*calpha)
	go to 999

C **********************************************************************
C ISCAT= 1 : Isotropic
20	p = 1.0d0
	go to 999

C **********************************************************************
C ISCAT= 2 : Henyey-Greenstein
30	f1 = cons(1)
	f2 = 1.0d0-f1
	hg11 = 1.0d0-cons(2)*cons(2)
	hg12 = 2.0d0-hg11
	hg21 = 1.0d0-cons(3)*cons(3)
	hg22 = 2.0d0-hg21
C        print*,hg11,hg12,hg21,hg22,calpha
C        print*,cons(1),cons(2),cons(3)
	p = f1 * hg11 / dsqrt( hg12 - 2.0d0*cons(2)*calpha )**3 +
     1   f2 * hg21 / dsqrt( hg22 - 2.0d0*cons(3)*calpha )**3
	go to 999

C **********************************************************************
C ISCAT= 3 : Legendre polynomial
40	IF (icons.GT.maxscatpar) then
          PRINT*,' PHASE1.f :: Error. NLEG is too large.'
          PRINT*,' PHASE1.f :: Stopping program.'
          STOP
        ENDIF
cc40	if (icons.gt.maxscatpar) call abend(' PHASE1:  NLEG TOO LARGE')
	cons(1) = 0.25d0/pi
	p = 0.d0
	xf = calpha*calpha
	do k=1,(icons+1)/2	! polynomials of order 2k-2 & 2k-1
	  n0 = 2*k-1
	  n1 = 2*k
	  x0 = 1.d0
	  x1 = calpha
	  pa0 = 0.d0
	  pa1 = 0.d0
	  do i = 1,k		! powers 2i-2 & 2i-1
	    pa0 = pa0 + x0*dfloat(jcoef(i,n0))/dfloat(jdiv(k))
	    pa1 = pa1 + x1*dfloat(jcoef(i,n1))/dfloat(jdiv(k))
	    x0 = x0*xf
	    x1 = x1*xf
	  enddo
	  p = p+cons(n0)*pa0+cons(n1)*pa1
	enddo
C This function is assumed properly normalized [by cons(1)], so:
	return

C **********************************************************************
C ISCAT= 4 : External phase file
50      if(vwave.ne.oldvwave)then 
         call interp_phase(vwave,ncont)
         oldvwave=vwave
        end if

        do i=1,nphas
         pfunc(i)=phas(icont,3,nphas - i + 1)
         xmu(i)=cthet(nphas - i + 1)
        end do

    

        calph1=sngl(calpha)
        call verint(xmu,pfunc,nphas,p1,calph1)
        p=dble(p1)
	return

C **********************************************************************
C External Henyey-Greenstein file
60      call get_hg(vwave,calpha,ncont,icont,p)

        goto 999

C **********************************************************************
C Normalization factor 
999	p = p/(4.0d0*pi)
	return

	end
