************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C_TITLE:                SUBROUTINE CSPLINT
C
C_DESC:	SPLINT (for SPLINe INTerpolation) takes the given the arrays
C	XA(1,...,N), and YA(1,...,N), which tabulate a function (with 
C	the XA elements in order) and the given array Y2A, which is the
C	output of spline.f, and given a value of X, this routine returns a
C	cubic spline interpolated value Y.
C	The goal of the cubic spline interpolation is to get an
C	interpolation formula that is smooth in the first derivative,
C	and continuous in the second derivative, both within the interval
C	and the boundaries. The solution is ...
C
C		Y = Ay(j) + By(j+1) + Cy"(j) + Dy"(j+1)
C
C                   x(j+1) - x                          x - x(j)
C	where A = ---------------,	B = 1 - A = ---------------,
C                  x(j+1) - x(j)                     x(j+1) - x(j)
C
C                 1
C             C = - (A^3 - A)(x(j+1) - x(j))^2,
C                 6
C
C                 1
C             D = - (B^3 - B)(x(j+1) - x(j))^2
C                 6
C
C_HIST:	???????	???	ORIGINAL VERSION
C	23/11/07 PGJI   Updated to double precision to prevent return
C			of unassigned variables
C	1/11/11  PGJI	Changed name and moved to rtm_util
C-----------------------------------------------------------------------

	SUBROUTINE CSPLINT(XA,YA,Y2A,N,X,Y)
        IMPLICIT NONE
        INTEGER N,K,KLO,KHI
	REAL XA(N),YA(N),Y2A(N),X,Y
        DOUBLE PRECISION Y1,A,B,H
C-----------------------------------------------------------------------
C
C	Bracket KLO and KHI around the input value X
C
C-----------------------------------------------------------------------

	KLO=1
	KHI=N
1	IF (KHI-KLO.GT.1) THEN
		K=(KHI+KLO)/2
		IF(XA(K).GT.X)THEN
			KHI=K
		ELSE
			KLO=K
		ENDIF
		GOTO 1
	ENDIF

C-----------------------------------------------------------------------
C
C	Determine whether the XA elements are distinct
C
C-----------------------------------------------------------------------

	H=DBLE(XA(KHI)-XA(KLO))
	IF (H.EQ.0.) THEN
		WRITE(*,*)'CSPLINT: Bad XA input.'
		WRITE(*,*)N,X
		WRITE(*,*)XA(1),XA(N)
		STOP
	ENDIF

C-----------------------------------------------------------------------
C
C	Evaluate the cubic-spline polynomial
C
C-----------------------------------------------------------------------

	A=DBLE(XA(KHI)-X)/H
	B=DBLE(X-XA(KLO))/H
	Y1=A*DBLE(YA(KLO))+B*DBLE(YA(KHI))+
     1  ((A**3-A)*DBLE(Y2A(KLO))+(B**3-B)*DBLE(Y2A(KHI)))*(H**2)/6.
        
        IF(ABS(Y1).GT.1.00E-37)THEN
         Y=SNGL(Y1)
        ELSE
         Y=1.00E-37
        ENDIF
C        print*,'csplint',Y1,Y
C-----------------------------------------------------------------------
C
C       Return and end.
C
C-----------------------------------------------------------------------

	RETURN

	END

************************************************************************
************************************************************************
