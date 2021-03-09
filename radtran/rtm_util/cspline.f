************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C_TITLE:		SUBROUTINE CSPLINE
C
C_DESC:	Given arrays x(1,...,N) and y(1,...,N) containing a tabulated
C	function, and given values for YP1 and YPN for the first
C	derivative of the interpolating function at points 1 and N,
C	respectively, this routine returns an array Y2(1,...,N) that
C	contains the second derivatives of the interpolating function at
C	the tabulated points x(1,...N). If YP1 and/or YPN are equal to
C	1e30 or larger, the routine is signaled to set the corresponding
C	boundary condition for a natural spline, with zero derivative on
C	that boundary.
C	This program need only be called once.
C
C_HIST:	???????	???	ORIGINAL VERSION
C	13.2.90	SBC	Error trap "IF(N.GT.NMAX)THEN" added.
C       1.2.11  PGJI    Changed name and set dimension of U from arrdef.f 
C-----------------------------------------------------------------------

	SUBROUTINE CSPLINE(X,Y,N,YP1,YPN,Y2)

	IMPLICIT NONE

C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)
        INCLUDE '../includes/arrdef.f'

	INTEGER		I,K,N
	REAL		YP1,YPN
	REAL		P,QN,UN,SIG

	REAL		X(N),Y(N),Y2(N),U(MAXBIN)
        integer idiag,iquiet
        common/diagnostic/idiag,iquiet

C-----------------------------------------------------------------------

	IF(N.GT.MAXBIN)THEN
	 WRITE(*,*)'Error in CSPLINE - recompile with larger MAXBIN'
	 WRITE(*,*)' N= ', n, ' maxbin= ', maxbin
	 STOP
	ENDIF

C-----------------------------------------------------------------------
C
C	Setting the lower boundary condition to be either equal to zero,
C	the so-called "natural cubic spline", or else to have a specified
C	first derivative.
C
C-----------------------------------------------------------------------

	IF(YP1.GT..99E30)THEN
		Y2(1)=0.
		U(1)=0.
	ELSE
		Y2(1)=-0.5
		U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
	ENDIF
 
C-----------------------------------------------------------------------
C
C	The following is the decomposition loop of the tridiagonal
C	algorithm. Y2 and U are used for temporary storage of the
C	decomposed factors.
C
C-----------------------------------------------------------------------


	DO 11 I=2,N-1
		SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
		P=SIG*Y2(I-1)+2.
		Y2(I)=(SIG-1.)/P

		IF(X(I-1).EQ.X(I))THEN
			U(I)= U(I-1)
			GOTO 11
		ELSEIF(X(I).EQ.X(I+1))THEN
			U(I)= U(I-1)
			GOTO 11
		ENDIF

		U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1		/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P

11    CONTINUE

C-----------------------------------------------------------------------
C
C	Setting the upper boundary condition to be either equal to zero,
C	the so-called "natural cubic spline", or else to have a specified
C	first derivative.
C
C-----------------------------------------------------------------------

	IF (YPN.GT..99E30) THEN
		QN=0.
		UN=0.
	ELSE
		QN=0.5
		UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
	ENDIF
	Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)

C-----------------------------------------------------------------------
C
C	The following is the backsubstituition loop of the tridiagonal
C	algorithm.
C
C-----------------------------------------------------------------------

	DO 12 K=N-1,1,-1
		Y2(K)=Y2(K)*Y2(K+1)+U(K)
12	CONTINUE

        do i=1,n
         if(isnan(y2(i)).and.idiag.gt.0)then
             print*,'Error in cspline.f - NAN returned.',I
         endif
        enddo
C-----------------------------------------------------------------------
C
C       Return and end.
C
C-----------------------------------------------------------------------

	RETURN

	END

************************************************************************
************************************************************************
