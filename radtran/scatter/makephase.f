C-----------------------------------------------------------------------
C
C				PROGRAM MAKEPHASE
C
C	Calculates phase functions of spherical particles according to 
C	Mie theory and outputs these to file.
C	Mie scattering computations use subroutine MIESCAT (& DMIE).
C
C	3 Mie size distributions are allowed:
C
C         (1)  Standard Gamma distribution
C              N(r) = const*r**ALPHA * exp( -r/(AA*BB))
C              Input constants are AA,BB ( ALPHA = (1-3*BB)/BB) 
C
C         (2)  Log-Normal distribution:
C              N(r) = 1/(sqrt(2*pi)*r*SIG)*exp(-((log(r)-log(R0))/SIG)**2/2)
C              Input constants are R0,SIG
C       
C         (3)  MCS Gamma distribution
C              N(r)=const * r**a * exp(-b*r**c)
C
C	In addition the following options are available:
C 	
C	  (4)  Mie phase function for a single particle size
C	  (5)  Isotropic scattering
C	  (6)  Henyey-Greenstein scattering
C	  (7)  Dipole scattering
C
C	The size integration parameters are RMIN, RMAX, DELTA(R)
C	(coded as RS(3)). If RMAX<RMIN, subroutine MIESCAT will 
C	terminate the integration at the point where all modes have 
C	N*QSCAT < 10**-6 * (N*QSCAT)max (this option is currently 
C	forced).
C	
C	The CPU time depends strongly on the series summation in subr. 
C	DMIE, which is terminated when the Mie coefficients become 
C	< 1.e-14, a circumstance that is not simple to predict.  However, 
C	some guidance may be obtained from the upper limit on the summation 
C	that is hard- coded into the subroutine:  this limit is constant 
C	(135) - hence the CPU time should be roughly proportional to the 
C	number of integration steps NR = (RMAX-RMIN)/DR -- as long as 
C	XMAX*M <= 135, where:
C
C		XMAX = 2*pi*RMAX/LAMBDA,
C		M = sqrt(RR*RR+RI*RI), and
C		(RR,RI) is the complex index of refraction.
C
C	However, if XMAX*M > 135, then the limit increases with RMAX, 
C	hence the CPU time may be expected to approach NR**2 for 
C	increasing RMAX.
C
C
C	PGJI	12/12/94Based on an original version by Kamp and Collard.
C	ALW	10/1/96	Further modified for use by NIMSRAD.
C       CAN	23/5/96	Took Andy's version and added Titan tholins
C       CAN	19/7/96	Added in Methane Ice spectral data
C	PGJI	16/2/12 Overhauled and extended to fit HG functions (i.e. 
C			to mimic the now defunct Ozonetable)
C-----------------------------------------------------------------------

	program makephase

	implicit none
	integer	max_theta, max_mode, max_wave,inorm
	parameter (max_theta = 100)
	parameter (max_mode = 10)
	parameter (max_wave = 1000)

	integer	nmode, iwave, nwave, ntheta, nphase, I, iscat, iref, 
     1		irec, J, siref(max_mode), siscat(max_mode), wc, K
	real	pi, wmin, wmax, wstep, minlam, lambda, parm(3), 
     1		rs(3), refind(2), w,  theta(max_theta), scat, 
     2		ext, phase(max_theta), omega, calpha,  x, qext, gsec,
     3		xsec(max_mode,max_wave,2), qsca, qabs, sum, du, dsum, 
     4		sparm(max_mode,3), srs(max_mode,3),
     5		scalef,srefind(max_mode,2),f,g1,g2,rms,
     6          phase1(max_theta)
	double precision f1, f2, hg11, hg12, hg21, hg22
	complex	nc

	character*1 	ans
	character*10	wavetype(2)
	character*100 	outfile,hgfile
	character*512 	buffer
      

        integer         ilist, nspec, idspec
        real          table(5000,3)
        common /store/ ilist,nspec,idspec,table


	pi = 3.141592654
	wavetype(1) = 'wavelength'
	wavetype(2) = 'wavenumber'

C-----------------------------------------------------------------------
C
C	Get input. 
C
C-----------------------------------------------------------------------

	write (*,*) '           '
	write (*,*) '	 Welcome to Makephase'
	write (*,*) '           '

	write (*,*) '           '
10	call prompt('Enter number of aerosol modes: ')
	read  (*,*) nmode
	if (nmode.gt.max_mode) then
		write (*,*) ' Maximum number of modes = ', max_mode
		goto 10
	endif

20	write (*,'('' Select (1) wavelength (microns) or (2)
     1 wavenumbers (cm-1): '',$)')
	read (*,*) iwave
	if ((iwave.ne.1).and.(iwave.ne.2)) goto 20

30	write (*,'('' Enter range, (start, end, delta): '',$)')
	read (*,*) wmin, wmax, wstep
	if (wmin.ge.wmax) then
		write (*,*) ' Incompatible min/max values.'
		goto 30
	endif

        write (*,'('' Enter name of .xsc file to be created: '',$)')
        read (*,'(a60)') outfile
	call file (outfile, outfile,'xsc')
        open (12, file=outfile, status='unknown')

        call prompt('Renormalise Phase Function (Y/N)? : ')
        read(5,1)ans
1       format(a)

        inorm=0
        if(ans.eq.'y'.or.ans.eq.'Y')inorm=1

C-----------------------------------------------------------------------
C
C	Initialise some stuff. Minlam is the minimum wavelength of the 
C	range, and is passed to GET_INPUT for potential use in 
C	selecting the size integration range.
C
C-----------------------------------------------------------------------

	nwave = 1 + int((wmax-wmin)/wstep)
	if (nwave.gt.max_wave) then
		write (*,*) ' Maximum number of steps exceeded'
		write (*,*) ' Maximum = ', max_wave, ' Current = ',nwave
		stop
	endif

	call get_theta (theta, max_theta, ntheta, nphase)

	if (iwave.eq.1) then
		minlam = wmin
	else
		minlam = 1.e4/wmax	! From wavenumbers
	endif

C-----------------------------------------------------------------------
C
C	Begin loop over number of modes. Get inputs per mode. These are
C	stored so that all inputs can be given near the beginning of 
C	the program run
C
C-----------------------------------------------------------------------

	do I = 1, nmode
		write (*,*) '          '
		write (*,*) ' Beginning inputs for mode', I
		write (*,*)
		
		call get_input(minlam, iscat, iref, parm, rs, scalef)
		siscat(I) = iscat
		siref(I) = iref
		do J = 1, 3
			sparm(I,J) = parm(J)
			srs(I,J) = rs(J)
		enddo

		if (iref.eq.1) then
			print*,'Give real and imaginary parts'
     			call prompt('of refractive index: ')
			read (*,*)  refind(1), refind(2)
			srefind(I,1) = refind(1)
			srefind(I,2) = refind(2)
		endif
	enddo

C-----------------------------------------------------------------------
C
C	Begin calculation over number of modes. Get relevant inputs. 
C	Write mode dependent headers to the PHASEN files (N is replaced 
C	by the number of the mode under consideration).
C
C-----------------------------------------------------------------------

	do I = 1, nmode
 		ilist=0
		iscat = siscat(I)
		iref = siref(I)
		do J = 1, 3
			parm(J) = sparm(I,J)
			rs(J) = srs(I,J)
		enddo

                print*,iscat,iref,parm,rs

		if (iref.eq.1) then
			refind(1) = srefind(I,1)
			refind(2) = srefind(I,2)
		endif

		outfile = 'PHASEN.DAT'
		outfile(6:6) = char(I+48)
		open (13,file=outfile, status='unknown', 
     1                 access='direct',recl=512, form='formatted')
		irec = 1
		write (buffer, 1010) wavetype(iwave), wmin, wmax, wstep,
     1				nwave, nphase
		write (13, 1000, rec=irec) buffer
		irec = irec + 1

		if (iscat.eq.1) then
			write (buffer, 1020) (parm(J), J = 1, 3),
     1				(rs(J), J = 1, 3)
		elseif (iscat.eq.2) then
			write (buffer, 1030) (parm(J), J = 1, 2),
     1				(rs(J), J = 1, 3)
		elseif (iscat.eq.3) then
			write (buffer, 1035) (parm(J), J = 1, 3),
     1				(rs(J), J = 1, 3)
		elseif (iscat.eq.4) then
			write (buffer, 1040) parm(1)
		elseif (iscat.eq.5) then
			write (buffer, 1050) 
		elseif (iscat.eq.6) then
			write (buffer, 1060) (parm(J), J = 1, 3)
		elseif (iscat.eq.7) then
			write (buffer, 1070) parm(1)
		endif

		write(13,1000,rec=irec) buffer
		irec = irec + 1

		write (buffer, 1080) (theta(J), J = 1, nphase)
		write(13,1000,rec=irec) buffer
		irec = irec + 1


		hgfile = 'hgphaseN.DAT'
		hgfile(8:8) = char(I+48)
		open (14,file=hgfile, status='unknown')
                

C-----------------------------------------------------------------------
C
C	Begin loop over wavenumbers or wavelengths.
C
C-----------------------------------------------------------------------

		do J = 1, nwave
			w = wmin + (J-1) * wstep
			if (iwave.eq.1) then
				lambda = w
			else
				lambda = 1.e4/w
			endif

			if (iref.eq.2) then
				write (*,'('' Give real and imaginary''
     1					'' parts of refractive'',
     2					'' index: '',$)')
				read (*,*)  refind(1), refind(2)
			elseif (iref.eq.3) then
				call h2o_refind(lambda, refind(1),
     1				       refind(2))
			elseif (iref.eq.4) then
			        call nh3_refind(lambda, refind(1),
     1				       refind(2))
			elseif (iref.eq.5) then
			        call tholin_refind(lambda, refind(1), 
     1                                 refind(2), scalef)
			elseif (iref.eq.6) then
			        call methane_refind(lambda, refind(1), 
     1                                 refind(2))
			elseif (iref.eq.7) then
			        call nh4sh_refind(lambda, refind(1), 
     1                                 refind(2))
			elseif (iref.eq.8) then
			        call hydra_refind(lambda, refind(1), 
     1                                 refind(2))
			elseif (iref.eq.9) then
			        call table_refind(lambda, refind(1), 
     1                                 refind(2))
			endif
	
			print*, 'Refractive index: ', refind(1),
     1				 refind(2)

C-----------------------------------------------------------------------
C
C	Do Mie scattering first. Omega is the single scattering albedo.
C
C-----------------------------------------------------------------------

			if (iscat.lt.5) then
				call miescat(lambda, iscat,parm, 1, rs,
     1					refind, theta, 
     2					ntheta, scat, ext, phase, 
     3				        nphase)
			omega = scat/ext

			end if


C-----------------------------------------------------------------------
C
C	Do isotropic scattering.
C
C-----------------------------------------------------------------------

			if (iscat.eq.5) then	
				do K = 1, nphase
					phase(K) = 1.0
         			end do
			endif
			if ((iscat.eq.5).or.(iscat.eq.6)) then
				if (iwave.eq.1) then
					write (*,* ) 
     1						' Wavelength is: ', w
				else
					write (*,* ) 
     1						' Wavenumber is: ', w
				endif
  				write (*,'('' Give extinction'',
     1	'' cross-section (cm2) and single scattering albedo: '',$)')
				read (*,*) ext, omega
			endif

C-----------------------------------------------------------------------
C
C	Do Henyey-Greenstein scattering.
C
C-----------------------------------------------------------------------

			if (iscat.eq.6) then
				f1 = parm(1)
				f2 = 1.0d0-f1
				hg11 = 1.0d0-parm(2)*parm(2)
				hg12 = 2.0d0-hg11
				hg21 = 1.0d0-parm(3)*parm(3)
				hg22 = 2.0d0-hg21
				do K = 1, nphase	
					calpha=cos(theta(K)*pi/180)
					phase(K) = sngl(f1 * hg11 / 
     1	dsqrt( hg12 - 2.0d0*parm(2)*calpha )**3 + f2 * hg21 / 
     2	dsqrt( hg22 - 2.0d0*parm(3)*calpha )**3)
				end do
			endif

C-----------------------------------------------------------------------
C
C	Do dipole scattering.
C
C-----------------------------------------------------------------------

			if (iscat.eq.7) then	
				nc=cmplx(refind(1),-refind(2))
				x = 2 * pi * parm(1)/lambda
				qsca = (8./3.)*(x**4) * abs((nc**2-1)/
     1					(nc**2 + 2))
				qabs = - 4 * x * aimag((nc**2 - 1)/
     1					(nc**2 + 2))
				qext = qsca + qabs
				omega = qsca/qext
				gsec = pi * (parm(1) * 1.e-4)**2
				ext = qext * gsec
   
				do K = 1, nphase
					calpha = cos(theta(K)*pi/180)
 					phase(K) = 0.75 * (1.0 + 
     1						calpha*calpha)
				end do
			end if 
  
C-----------------------------------------------------------------------
C
C	Write output data at each wavenumber or wavelength, then close 
C	PHASEN file. Also perform renormalization to ensure that the
C       intgeral of the phase function from cos(theta) = -1 to +1 is equal
C       to 1/2pi
C
C-----------------------------------------------------------------------

	  	 	do K = 1, nphase
				phase1(K)=phase(k)
			 	phase(K) = 0.25 * phase(K)/pi
		        enddo


                        sum=0
                        do k=1,nphase-1
                         du = cos(theta(k)*pi/180.0) - 
     &                        cos(theta(k+1)*pi/180.0)
                         sum=sum+0.5*(phase(k)+phase(k+1))*du
                        end do
                        sum=sum*2*pi
                        print*,'Makephase: 2*pi*Int(pdu) = ',sum
                        dsum = abs(sum - 1.0)
                        if(dsum.gt.0.05)then
                         print*,'Makephase: WARNING, > 5% off unity'       
                        end if
                        if(inorm.eq.1)then
                         print*,'Renormalising : '
                         do k=1,nphase
                          phase(k)=phase(k)*1.0/sum
                          phase1(k)=phase1(k)*1.0/sum
                         end do
                        endif

                        if(phase(1).ge.10000.0)then
 			 write (buffer,1095,ERR=666) 
     &			   w, ext, omega,(phase(K),K=1, nphase)
			else
			 write (buffer,1090,ERR=666) 
     &			   w, ext, omega,(phase(K),K=1, nphase)
			endif
        		write(13,1000,rec=irec)buffer

			goto 667
 666			print *, 'WARNING: error occured during write',
     &                       ' to buffer. Probably the forward'
			print *, 'scattering amplitude exceeded the',
     &                    ' FORMAT limits due to large particle size.'

 667			irec=irec+1

C			Subfithgm is expecting phase functions normalised
C			to 4pi, rather than 1, so pass phase1 instead of
C			phase.
			call subfithgm(nphase,theta,phase1,f,g1,g2,rms)

                        write(14,*)w,f,g1,g2

			xsec(I,J,1) = ext
			xsec(I,J,2) = omega
		enddo

 		close(13)
 		close(14)
	enddo

C-----------------------------------------------------------------------
C
C	Now produce .XSC files
C
C-----------------------------------------------------------------------

	write (12,*) nmode        
	do I = 1, nwave
		w = wmin + (I-1) * wstep
		write(12,*) w,(xsec(J,I,1), J=1, nmode)
		write(12,*)(xsec(J,I,2), J=1, nmode)
	enddo

        close(12)

C-----------------------------------------------------------------------
C
C	Format statements and end.
C
C-----------------------------------------------------------------------

1000	format (a512)
1010	format (1x, a10, 2(2x, f8.2), 2x, f8.4, 2(2x, i4))
1020	format (' Mie, standard gamma, A = ',0pf7.2, ' B = ',0pf7.2,
     1		' alpha = ',0pf7.2, ' rmin = ', 0pf7.3, ' rmax = ',
     2		0pf7.2, ' rstep = ', 0pf7.3)
1035	format (' Mie, MCS gamma, A = ',0pf7.2, ' B = ',0pf7.2,
     1		' C = ',0pf7.2, ' rmin = ', 0pf7.3, ' rmax = ',
     2		0pf7.2, ' rstep = ', 0pf7.3)
1030	format (' Mie, log normal ', ' R0 = ',0pf7.2, ' Sig = ',0pf7.2,
     1		' rmin = ', 0pf7.3, ' rmax = ',0pf7.2, ' rstep = ', 
     2		0pf7.3)	
1040	format (' Mie, single particle ', ' R0 = ',0pf7.2)	
1050	format (' Isotropic scattering')
1060	format (' Henyey-Greenstein, f = ', 0pf7.2, ' g1 = ',0pf7.2,
     1		' g2 = ',0pf7.2)
1070	format (' Dipole, R0 = ', 0pf7.2)
1080	format (100(1x,f8.3))
1090    format (1x,f8.2,1x,e12.5,1X,e12.5,50(1x,f10.5))
1095    format (1x,f8.2,1x,e12.5,1X,e12.5,50(1x,f10.4))

	end


*************************************************************************
*************************************************************************
C-----------------------------------------------------------------------
C
C		SUBROUTINE GET_INPUT
C
C	Delivers input variables for Makephase for each mode.
C
C-----------------------------------------------------------------------

	subroutine get_input(minlam, iscat, iref, parm, rs, scalef)

	implicit none
	integer	iscat,  isize, iref, I
	real	minlam, parm(3), rs(3), scalef

10      print*,'Calculation types are : '
        print*,'(1) Mie scattering, standard gamma distribution'
        print*,'(2) Mie scattering, log-normal distribution'
        print*,'(3) MCS Modified Gamma distribution'
        print*,'(4) Mie scattering,  single particle size'
        print*,'(5) Isotropic scattering'
        print*,'(6) Henyey-Greenstein scattering'
        print*,'(7) Dipole scattering'

        call prompt('Select one : ')
	read (*,*) iscat
	if ((iscat.lt.1).or.(iscat.gt.7)) goto 10

	write (*,*) '           '

C-----------------------------------------------------------------------
C
C	Get size distribution parameters.
C
C-----------------------------------------------------------------------

	do I = 1, 3
		parm(I) = 0.
	enddo

	if (iscat.eq.1) then
		print*,'Give A and B (effectively, A is the '
                print*,'mean radius of the distribution in'
                call prompt('microns, and B is the variance) : ')
		read(*,*) parm(1), parm(2)
		parm(3) = (1. - 3 * parm(2))/parm(2)
		write (*,*) '           '
	elseif (iscat.eq.2) then
		print*,'Give R0 and SIG (effectively, R0 is the'
                print*,'mean radius of the distribution in'
                print*,'microns, and SIG is the mean of '
                call prompt('log(r)/log(R0)) : ')
		read(*,*) parm(1), parm(2)
		write (*,*) '           '
	elseif (iscat.eq.3) then
		call prompt('Give A, B and C : ')
		read(*,*) parm(1), parm(2), parm(3)
		write (*,*) '           '
	elseif ((iscat.eq.4).or.(iscat.eq.7)) then
		print*,'Give R0 (R0 is the particle size in '
     		call prompt('microns) : ')
		read (*,*) parm(1)
	endif	

C-----------------------------------------------------------------------
C
C	Get size integration parameters (if required). If RMAX < RMIN, 
C	subroutine MIESCAT will terminate the integration at the
C	point where all modes have N*QSCAT < 10**-6 * (N*QSCAT)max.
C	For one size Mie particles, the size distribution parameters are 
C	set to those of the particle.
C
C-----------------------------------------------------------------------

	isize = 0
	do I = 1, 3
		rs(I) = 0.
	enddo

	if (iscat.lt.4) then
20		print*,'Either specify size integration '
             call prompt('parameters (1) or use default (2) : ')
		read (*,*) isize
		if ((isize.ne.1).and.(isize.ne.2)) goto 20
	elseif (iscat.eq.4) then
		do I = 1, 3
			rs(I) = parm(1)
		enddo
	endif

	if (isize.eq.1) then
		call prompt('Give min, max and stepsize : ')
		read (*,*)  rs(1), rs(2), rs(3)
C		Don't know why this line was commented out.
C		rs(2) = 0.

		if ((rs(2)-rs(1))/rs(3).gt.10000) then
			write (*,*) ' Warning: Number of size steps',
     1				' exceeds 10000'

		endif

	elseif (isize.eq.2) then
		rs(1) = 0.015 * minlam
		rs(2) = 0.
		rs(3) = rs(1)
		write (*,*) 'Parameters used (min, max, and step):'
		write (*,*) '           ', rs

			write (*,*) 'Size integration will be',
     1			' terminated at n*Qscat=1.e-6*max(n*Qscat)'

	endif

C-----------------------------------------------------------------------
C
C	Get Henyey-Greenstein constants (if required) or get
C	refractive indices.
C
C-----------------------------------------------------------------------

	write (*,*) '          '
	if (iscat.eq.6) then
		print*,'Give 3 Henyey_Greenstein constants'
                call prompt('(f, g1, g2) : ')
		read (*,*) parm(1), parm(2), parm(3)
		iref = 0
	elseif (iscat.eq.5) then
		iref = 0
	else
30		print*,'Refractive index options are: '
                print*,'(1) Constant value over range'
                print*,'(2) Manual input at each calculation point'
     		print*,'(3) Water lookup table'
     		print*,'(4) Ammonia look up table'
                print*,'(5) Titan tholins'
            	print*,'(6) Methane look-up table'
     		print*,'(7) NH4SH look-up table'
     		print*,'(8) Hydrazine look-up table'
     		print*,'(9) read from external RI table'
     		call prompt('Select one : ')
		read (*,*) iref
		if ((iref.lt.1).and.(iref.gt.9)) goto 30
	endif

C allow optional scaling of tholin refractive index to match Titan
C haze:

	if (iref.eq.5) then
	   write (*,*) 'Enter scale factor for tholin imaginary',
     &                 ' refractive index:'
	   read(*,*) scalef
	endif

	return

	end


