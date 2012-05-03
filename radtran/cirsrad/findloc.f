************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			   SUBROUTINE FINDLOC
C
C	For a given monotonic array and 2 values, finds the indices of 
C	the array values contained by the passed variables. Uses 
C	numerical recipes subroutines.
C
C-----------------------------------------------------------------------
	
	subroutine findloc (x, nx, y1, y2, I1, I2)

	implicit none
	integer		nx, I1, I2, Ibot, Itop, Itmp, I, I1f, I2f
	real		x(nx), y1, y2, xbot, xtop

	if (x(1).lt.x(nx)) then
		xbot = x(1)
		xtop = x(nx)
		Ibot = 1
		Itop = nx
	else
		xbot = x(nx)
		xtop = x(1)
		Ibot = nx
		Itop = 1
	endif

C-----------------------------------------------------------------------
C
C	Check ranges are valid. Then for each y, do the following.
C	If one end of range is outside of array, assign array end values, 
C	else LOCATE returns Itmp such that y is between x(Itmp) and 
C	x(Itmp + 1). Check for an exact match, and flag if this happens
C
C-----------------------------------------------------------------------

	if (((y1.lt.xbot).and.(y2.lt.xbot)).or.((y1.gt.xtop).and.
     1		(y2.gt.xtop))) then
		I1 = 0
		I2 = 0
		goto 10
	endif

	if (y1.le.xbot) then
		I1 = Ibot
		I1f = 1
	elseif (y1.ge.xtop) then
		I1 = Itop
		I1f = 1
	else
		call locate (x, nx, y1, Itmp)
		I1 = Itmp
		I1f = 0
		do I = Itmp, Itmp+1
			if (abs(y1-x(I)).lt.1.e-5) then
				I1 = I
				I1f = 1
			endif
		enddo
	endif

	if (y2.le.xbot) then
		I2 = Ibot
		I2f = 1
	elseif (y2.ge.xtop) then
		I2 = Itop
		I2f = 1
	else
		call locate (x, nx, y2, Itmp)
		I2 = Itmp
		I2f = 0
		do I = Itmp, Itmp+1
			if (abs(y2-x(I)).lt.1.e-5) then
				I2 = I
				I2f = 1
			endif
		enddo
	endif

C-----------------------------------------------------------------------
C
C	Make sure that the index for the array value above the bottom
C	of the range is returned. If there is an exact match, leave it
C	unaltered.		
C
C-----------------------------------------------------------------------

	if (y1.lt.y2) then
		if (I1f.eq.0) I1 = I1 + 1
	else
		if (I2f.eq.0) I2 = I2 + 1
	endif
		
10	return
		 	
	end

************************************************************************
************************************************************************
