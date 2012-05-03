      SUBROUTINE CUBINT(XA1,YA1,N,X,Y)
C     $Id: cubint.f,v 1.1.1.1 2000-08-17 09:26:56 irwin Exp $
C     ************************************************************
C     Numerical recipes routine to do cubic spline interpolation
C     of a table
C
C     Pat Irwin		26/3/98
C
C     ************************************************************
      PARAMETER (NMAX=1000) 
      DIMENSION XA1(N),YA1(N)
      DIMENSION XA(NMAX),YA(NMAX),Y2(NMAX),Y2H(NMAX)
      DIMENSION XH(NMAX),YH(NMAX)
      LOGICAL RECALC
      COMMON /CHOLD/ XH,YH,Y2H

      IF(XA1(1).GT.XA1(2))THEN
       DO I=1,N
        XA(I) = XA1(1+N-I)
        YA(I) = YA1(1+N-I)
       ENDDO
      ELSE
       DO I=1,N
        XA(I) = XA1(I)
        YA(I) = YA1(I)
       ENDDO
      ENDIF

      RECALC=.FALSE.
      DO I=1,N
       IF(XA(I).NE.XH(I).OR.YA(I).NE.YH(I))RECALC=.TRUE.
      ENDDO

      IF(RECALC)THEN
C       print*,'cubint. Calculating new derivatives'
       YP1 = 0.0
       YPN = 0.0
       CALL CSPLINE(XA,YA,N,YP1,YPN,Y2)
       DO I=1,N
        XH(I) = XA(I)
        YH(I) = YA(I)
        Y2H(I) = Y2(I)
       ENDDO
      ELSE
C       print*,'Using stored derivatives'
       DO I=1,N
        Y2(I) = Y2H(I)
       ENDDO
      ENDIF

      CALL CSPLINT(XA,YA,Y2,N,X,Y)

      RETURN

      END
