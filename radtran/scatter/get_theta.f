*************************************************************************
*************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE GET_THETA
C
C	Calculates scattering angle array.
C
C	Theta		Scattering angles
C	Max_theta	Maximum size of theta array
C	ntheta		Size of theta subarray where theta .le. 90
C
C-----------------------------------------------------------------------

	subroutine get_theta (theta, max_theta, ntheta, nphase)

	implicit none
	integer	max_theta, ntheta, nphase, I
	real	theta(max_theta), dtheta

	ntheta = 1
	theta(1) = 0.
	dtheta = 1.0

10	if ((theta(ntheta)+dtheta).le.90.) then
		ntheta = ntheta + 1
		theta(ntheta) = theta(ntheta-1) + dtheta
		if (theta(ntheta).ge.5.) dtheta = 2.5
		if (theta(ntheta).ge.20.) dtheta = 5.0
		if (theta(ntheta).ge.40.) dtheta = 10.
		goto 10
	endif

	nphase = 2*ntheta - 1
	do I = ntheta+1, nphase
		theta(I) = 180.- theta(nphase+1-I)
	enddo

	return

	end
