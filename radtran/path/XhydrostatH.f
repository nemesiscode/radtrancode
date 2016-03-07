      SUBROUTINE XHYDROSTATH(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1 IDGAS,ISOGAS,H,P,T,VMR,SCALE)
C     ****************************************************************
C     Subroutine to rescale the heights of a H/P/T profile according to
C     the hydrostatic equation above and below the level where height=0.
C
C     Pat Irwin		27/8/06
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'

      INTEGER NPRO,IPLANET,JZERO,I,J,NVMR,AMFORM
      INTEGER ISOGAS(MAXGAS),IDGAS(MAXGAS)
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),LATITUDE,DELH,X,TEMP
      REAL RADIUS,G,MOLWT,R,SCALE(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL XVMR(MAXGAS),XMOLWT,CALCMOLWT,SH,ATDEPTH,ATDEPTH1
      REAL XDEPTH
      CHARACTER*8 PNAME


C     **************** Code **************

C     R no included in constdef.f (RGAS) and needs scaling for units here.
      R = RGAS*0.001

999   CONTINUE
      ATDEPTH = H(NPRO)-H(1)

C     First find level closest to zero altitude
      DELH=1000.0
      DO I=1,NPRO
       X = ABS(H(I))
       IF(X.LT.DELH)THEN
        DELH=X
        JZERO=I
       ENDIF

       DO J=1,NVMR
        XVMR(J)=VMR(I,J)
       ENDDO

       CALL NEWGRAV(IPLANET,LATITUDE,H(I),RADIUS,G,PNAME)
       IF(AMFORM.EQ.0)THEN
        XMOLWT=MOLWT
       ELSE
        XMOLWT=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
       ENDIF

       SCALE(I)=R*T(I)/(XMOLWT*G)

      ENDDO

      print*,'XhydrostatH: JZERO = ',JZERO

      IF(JZERO.GT.1.AND.JZERO.LT.NPRO)THEN
       H(JZERO)= 0.0
      ENDIF

      DO 301 I=JZERO+1,NPRO
           SH=0.5*(SCALE(I-1)+SCALE(I))
           H(I) = H(I-1) - SH*LOG(P(I)/P(I-1))
301   CONTINUE

      DO 311 I=JZERO-1,1,-1
           SH=0.5*(SCALE(I)+SCALE(I+1))
           H(I) = H(I+1) - SH*LOG(P(I)/P(I+1))
311   CONTINUE

      ATDEPTH1=H(NPRO)-H(1)
   
      XDEPTH = 100*ABS((ATDEPTH1-ATDEPTH)/ATDEPTH)
      IF(XDEPTH.GT.1.0)THEN
        GOTO 999
      ENDIF

      RETURN
      END
