      SUBROUTINE DUSTPATH(PVEC,AVEC,NPRO,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,CONT)
C     $Id:
C     ***************************************************************
C     Subroutine to calculate equivalent CG path between two points
C     in a  spherically symmetric atmosphere. Code is based on the
C     RADTRAN layer.f routine and uses Simpson's rule integration.
C
C     Input variables:
C	PVEC(3)		REAL	Starting position. Convention is that
C				z-azis is the zenith at the tangent point
C				of a limb path.
C	AVEC(3)		REAL	Vector of new path.
C	NPRO		INTEGER	Number of points in atm profile.
C	NCONT		INTEGER	Number of dust types 
C	P(MAXPRO)		REAL	Pressure (prf)
C	T(MAXPRO)		REAL	Temperature (prf)
C	H(MAXPRO)		REAL	Heights (prf)
C	DUST(MAXPRO,MAXCON)REAL	Dust profiles
C
C     Output variables:
C	CONT(MAXCON)	REAL	Path dust amounts
C
C     Pat Irwin		1/4/99
C 
C     ***************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'
      REAL PVEC(3),AVEC(3),TVEC(3)
      INTEGER NPRO,NINT,J,NCONT,I,K,IFL
      REAL DUDS,DNOW,MOLWT,F
      REAL TNOW,PNOW,CALCALT,RADIUS,HEIGHT
      PARAMETER (NINT=101)
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO)
      REAL S,DELS,W(NINT)
      REAL TOTAM
      REAL DUST(MAXPRO,MAXCON),CONT(MAXCON),XSEC

C     Initialise some variables
C     ------------------------------------------------------------------
      TOTAM = 0.0
      
      DO 31 I=1,MAXCON
       CONT(I)=0.0
31    CONTINUE


      DO 20 I=1,NINT
       W(I)=2.
       IF(I.EQ.2*(I/2))W(I)=4.
20    CONTINUE
      W(1)=1.
      W(NINT)=1.
C     ------------------------------------------------------------------

C     Check some numbers
C     ------------------------------------------------------------------
      IF(NCONT.GT.MAXCON)THEN
       PRINT*,'ODPATH: NCONT > MAXCON'
       PRINT*,NCONT,MAXCON
       STOP
      ENDIF

      IF(NPRO.GT.MAXPRO)THEN
       PRINT*,'ODPATH: NPRO > MAXPRO'
       PRINT*,NPRO,MAXPRO
       STOP
      ENDIF

      S = 0.0
      DO 10 I=1,3
       S = S + AVEC(I)**2
10    CONTINUE

      S = SQRT(S)

 
      DELS=S/FLOAT(NINT-1)

      DO 30 I=1,NINT
       DO 25 J=1,3
        TVEC(J) = PVEC(J) + FLOAT(I-1)*AVEC(J)/FLOAT(NINT-1)
25     CONTINUE
       HEIGHT = CALCALT(TVEC,RADIUS)
       F=-1.0
       IFL = 0
       DO K=1,NPRO-1
        IF(HEIGHT.GE.H(K).AND.HEIGHT.LT.H(K+1))THEN
         F = (HEIGHT - H(K))/(H(K+1)-H(K))
         IFL=K
        ENDIF
       ENDDO
       IF(IFL.EQ.0)THEN
        IF(HEIGHT.LT.H(1))THEN
         IFL=1
         F=0.0
        ENDIF
        IF(HEIGHT.GE.H(NPRO))THEN
         IFL=NPRO-1
         F=1.0
        ENDIF
       ENDIF


       PNOW = (1-F)*P(IFL) + F*P(IFL+1)
       TNOW = (1-F)*T(IFL) + F*T(IFL+1)


C      calculating the number of molecules per km per cm2
       DUDS = MODBOLTZ*PNOW/TNOW
       TOTAM = TOTAM + DUDS*W(I)

       DO 124 J=1,NCONT
         DNOW = (1-F)*DUST(IFL,J) + F*DUST(IFL+1,J)
         CONT(J)=CONT(J)+DNOW*DUDS*W(I)*MOLWT/AVAGAD
124    CONTINUE

30    CONTINUE

      IF(TOTAM.NE.0.0)THEN
       DO 145 J=1,NCONT
        CONT(J)=CONT(J)*DELS/3
145    CONTINUE
      ENDIF

      RETURN

      END
