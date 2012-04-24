************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				SUBROUTINE CLOSE_SCAT
C
C	Designed to close any opened scattering files before next 
C	iteration of retrieval algorithm.
C
C       Pat Irwin	4/9/11	Original version
C	Pat Irwin	29/2/12	Updated for Radtrans2.0
C-----------------------------------------------------------------------

	subroutine close_scat(ncont)

	implicit none

	include '../includes/arrdef.f'

	integer	ncont, nmu, isol, lowbc, liscat(maxcon), 
     1          lnorm(maxcon), 
     1		lncons(maxcon), imu0, imu, nf, I, irec1(maxcon),
     2		nphas, maxrec(maxcon), J, K
	real	dist, aphi, lcons(maxcon,maxscatpar), thet(maxphas), 
     1		cthet(maxphas), phas(maxcon,3,maxphas), 
     2		head(maxcon,3,3)

	common/scatter/nmu, isol, dist, lowbc, liscat, lnorm, 
     1		lncons, lcons, imu0, imu, aphi, nf
        common/phase/thet,cthet,phas,head,irec1,nphas,maxrec

C-----------------------------------------------------------------------
C
C	Close file based on test.
C
C-----------------------------------------------------------------------

	do I = 1, ncont
		if (liscat(I).gt.3) close (10+I)
	enddo

C-----------------------------------------------------------------------
C
C	Reinitialise certain common block variables so old values 
C	aren't carried over.
C
C-----------------------------------------------------------------------

	do I = 1, maxcon
		irec1(I) = 0
		do J = 1, 3
			do K = 1, 3
				head(I,J,K) = 0.
			enddo
		enddo
	enddo

C-----------------------------------------------------------------------
C
C	Return and end
C
C-----------------------------------------------------------------------
	

	return

	end		

************************************************************************
************************************************************************
