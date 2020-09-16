      SUBROUTINE CALC_PCOEFF(NX,Y,X,XMIN,XMAX,COEFF)
C     $Id: calc_pcoeff.f,v 1.2 2011-06-17 15:40:25 irwin Exp $
C     ***********************************************************************
C     Calculates the polynomial fit of order IORDER (IORDP1-1) to the set
C     of points Yn(Xn) in the range XMIN to XMAX. 
C     The arrays Y and X are examined to see if the X-range is fully covered
C     and then polynomial coefficients COEFF are calculated.
C     If there are more than IORDER+2 points in the range, the coefficients
C     are calculated by a least squares fit using the Numerical Recipes
C     routine LFIT. If there are less than these number of points, the
C     highest possible order polynomial is calculated using the exact
C     matrix solution.
C
C     Input Variables
C	NX	INTEGER	Array dimension of Y and X.
C	X(NX)	REAL	Values of X for which Y is defined
C	Y(NX)	REAL	Corresponding values of Y
C	XMIN	REAL	Minimum of X-range for fit
C	XMAX	REAL	Maximum of X-range for fit
C     Output Variable
C	COEFF(IORDP1)	REAL	Calculated polynomial coefficients.
C
C_History:	07Oct94	 PGJI	Original. Adapted from code in GENLBL
C 
C     ***********************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NMAX,IMIN,IMAX,L,NX,NPT,ISIG,MORD,J,K,I
      PARAMETER (NMAX=100,MORD=10)
      REAL X(NX),Y(NX),X1(NMAX),COEFF(IORDP1),XMIN,XMAX
      DOUBLE PRECISION AA(MORD,MORD),BB(MORD,MORD)

C     Variables used by fitting routine LFIT
      REAL SIG(NMAX),COVAR(MORD,MORD),CHISQ
      INTEGER LISTA(MORD)
C------------------------------------------------------------------------------


      IF(NX.LT.2)THEN
       WRITE(*,111)
111    FORMAT(' %CALC_PCOEFF: Need more points in y(x) (nx < 2)')
       STOP
      END IF

C     ----------------------------------------------------------------
C     Find lowest point:
      DO 113 L=1,NX
       IF(X(L).GT.XMIN)GOTO 114
113   CONTINUE
      WRITE(*,115)
115   FORMAT(' %Calc_pcoeff: Y(X) doesn"t cover computation range')
      STOP

114   IMIN=L-1
      IF(IMIN.LT.1)IMIN=1
C     ----------------------------------------------------------------
C     Find upper point:
      DO 117 L=NX,1,-1
       IF(X(L).LT.XMAX)GOTO 118
117   CONTINUE
      WRITE(*,115)
      STOP

118   IMAX=L+1
      IF(IMAX.GT.NX)IMAX=NX
C     ----------------------------------------------------------------
C     Check number of ordinates and see if they are in range.

C      print*,'imin,imax = ',imin,imax

      IF((IMAX-IMIN+1).GT.NMAX)THEN
       WRITE(*,312)
312    FORMAT(' %Calc_pcoeff: filter grid too fine - increase NMAX')
       STOP
      END IF

C     ----------------------------------------------------------------
C     Calculate the relative X-scale for polynomial fitting
      NPT=0
      DO 116 L=IMIN,IMAX
       NPT=NPT+1
       X1(NPT)=X(L)-XMIN
C       print*,'x1,y',x1(npt),y(l)
116   CONTINUE

C      print*,'npt = ',npt

C     ----------------------------------------------------------------
C     Calculate the fit
      IF(NPT.GE.IORDER+2)THEN

C       Here if a least squares fit must be made to points over the bin
        K=IORDP1
        DO 127 ISIG=1,NPT
         SIG(ISIG)=1.
127     CONTINUE
        DO 126 L=1,K
         LISTA(L)=L
126     CONTINUE

        CALL LFIT(X1,Y(IMIN),SIG,NPT,COEFF,K,LISTA,K,COVAR,K,CHISQ)

      ELSE

C       here if insufficient points to fit polynomial so simply computing 
C       a NPT-1 order polynomial
        DO 121 L=1,NPT
         AA(1,L)=1.
         DO 122 J=2,NPT
          AA(J,L)=AA(J-1,L)*X1(L)
122      CONTINUE
121     CONTINUE

        CALL DMATINV(AA,NPT,MORD,BB)

        DO 123 L=1,IORDP1
         COEFF(L)=0.
123     CONTINUE

        DO 124 L=1,NPT
         DO 125 J=1,NPT
          COEFF(L)=COEFF(L)+SNGL(BB(J,L))*Y(IMIN-1+J)
125      CONTINUE
124     CONTINUE

      END IF

      RETURN

      END

