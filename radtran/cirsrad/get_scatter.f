************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				SUBROUTINE GET_SCATTER
C
C-----------------------------------------------------------------------

	subroutine get_scatter(radfile, ncont)

	implicit none

	include '../includes/arrdef.f'

	integer	nmu, isol, lowbc, liscat(maxcon), lnorm(maxcon), 
     1		lncons(maxcon), nf, npoint1, nphas, 
     2		maxrec(maxcon), iunit, irec1(maxcon), I, J,
     3		nc, ncont
	real	dist, aphi, v0, v1, dv, thet(maxphas), cthet(maxphas),
     1		phas(maxcon,3,maxphas), lcons(maxcon,maxscatpar),
     2		head(maxcon,3,3), sol_ang, emiss_ang
	double precision	mu1(maxmu), wt1(maxmu), galb

	character*10	wavetype
	character*30	scatfile(maxcon)
	character*100	radfile
	character*512	buffer

C - Common blocks

	common/scatd/mu1, wt1, galb 
	common/scatter1/nmu, isol, dist, lowbc, liscat, lnorm, 
     1		lncons, lcons, sol_ang, emiss_ang, aphi, nf
	common/phase/thet,cthet,phas,head,irec1,nphas,maxrec

C-----------------------------------------------------------------------
C
C	Read scattering information file.
C
C-----------------------------------------------------------------------

	WRITE(*,*)'calling read_scatter1'
	call read_scatter1(radfile,nmu,mu1,wt1,isol,dist,lowbc, galb, 
     1		nc, liscat, lnorm, lncons, lcons, scatfile, sol_ang, 
     2		emiss_ang, aphi, nf) 
	WRITE(*,*)'read_scatter1 OK'
        do i=1,nmu
         print*,i,mu1(i),wt1(i)
        enddo
C-----------------------------------------------------------------------
C
C	Check some variable sizes.
C
C-----------------------------------------------------------------------

	if (nmu.gt.maxmu) then
		write (*,*) ' GET_SCATTER: Too many mu points'
		write (*,*) ' Nmu = ',nmu,' Maxmu = ',maxmu
		stop
	endif

	if (nc.ne.ncont) then
		write (*,*) ' GET_SCATTER: Scatter dust continua',
     1			' not equal to path  continua'
		write (*,*) ' Scatter # = ',nc,' Path # = ', ncont
		stop
	endif
	do I = 1, ncont
		if (lncons(I).gt.maxscatpar) then
			write (*,*) ' GET_SCATTER: Too many constants',
     1				' defined'
			write (*,*) ' I = ', I, ' lncons = ', lncons(I),
     1				 ' Maxscatpar = ', maxscatpar
			stop
		endif
	enddo

C-----------------------------------------------------------------------
C
C	If scatfiles are supplied, these must be opened now for direct 
C	access by scattering subroutines.
C
C-----------------------------------------------------------------------

	do I = 1, nc
		if (liscat(I).eq.4) then 
			iunit = 10 + I
			open (iunit,file=scatfile(I),status='old', 
     1				access='direct', recl=512, form=
     2				'formatted')
			read (iunit,1000,rec=1) buffer
			read(buffer,1010) wavetype, v0, v1, dv, 
     1				npoint1, nphas
			maxrec(I) = npoint1 + 3
			irec1(I) = 3

			if (nphas.gt.maxphas) then
				write (*,*) ' GET_SCATTER: Too many',
     1					' phase points'
				write (*,*) ' I = ',I,' Maxphas = ',
     1					maxphas,' Nphas = ', nphas
				stop
			endif

			if (I.eq.1) then
				read(iunit,1000,rec=3) buffer
				read (buffer,*) (thet(J), J = 1, nphas)
				do J = 1, nphas
				cthet(J) = cos(thet(J)*0.0174532)
				enddo
			endif
		endif
	enddo

C-----------------------------------------------------------------------
C
C     	Return and end
C
C-----------------------------------------------------------------------

1000	format (a)
1010    format (1x, a10, 2(2x, f8.2), 2x, f8.4, 2(2x, i4))

	return

	end
		
************************************************************************
************************************************************************
