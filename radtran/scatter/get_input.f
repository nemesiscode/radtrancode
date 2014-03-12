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
     		print*,'(9) H2SO4 look-up table'
     		print*,'(10) Mars dust look-up table'
     		print*,'(11) read from external RI table'
     		call prompt('Select one : ')
		read (*,*) iref
		if ((iref.lt.1).and.(iref.gt.11)) goto 30
	endif

C allow optional scaling of tholin refractive index to match Titan
C haze:

	if (iref.eq.5) then
	   write (*,*) 'Enter scale factor for tholin imaginary',
     &                 ' refractive index:'
	   read(*,*) scalef
	endif

	if (iref.eq.9) then
	   write (*,*) 'Enter %acid concentration (between 75',
     &                 ' and 96% :'
	   read(*,*) scalef
	endif

	return

	end


