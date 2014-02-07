      SUBROUTINE XHYDROSTATP(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1 IDGAS,ISOGAS,H,P,T,VMR,HTAN,PTAN,SCALE)
C     ****************************************************************
C     Subroutine to rescale the pressures of a H/P/T profile according to
C     the hydrostatic equation above and below a specified altitude
C     given the pressure at that altitude.
C
C     Pat Irwin		27/8/06
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'

      INTEGER NPRO,IPLANET,JTAN,I,J,AMFORM,NVMR
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),LATITUDE,TEMP,SH
      REAL RADIUS,G,MOLWT,R,SCALE(MAXPRO),HTAN,PTAN,DELH
      REAL VMR(MAXPRO,MAXGAS)
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS)
      REAL XVMR(MAXGAS),XMOLWT,CALCMOLWT
      CHARACTER*8 PNAME


C     **************** Code **************

C     R now included in constdef.f (RGAS) and needs scaling for units here.
      R = RGAS*0.001

C     First find level immediately below the reference altitude
      JTAN=1
      DO I=1,NPRO
       IF(H(I).LT.HTAN)THEN
        JTAN=I
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

      SH = 0.5*(SCALE(JTAN)+SCALE(JTAN+1))
      DELH = H(JTAN+1)-HTAN
      P(JTAN+1)=PTAN*EXP(-DELH/SH)
      DELH = H(JTAN)-HTAN
      P(JTAN)=PTAN*EXP(-DELH/SH)

      DO 301 I=JTAN+2,NPRO
           SH = 0.5*(SCALE(I-1) + SCALE(I))
           DELH = H(I)-H(I-1)
           P(I)=P(I-1)*EXP(-DELH/SH)
301   CONTINUE

      DO 311 I=JTAN-1,1,-1
           SH = 0.5*(SCALE(I) + SCALE(I+1))
           DELH = H(I)-H(I+1)
           P(I)=P(I+1)*EXP(-DELH/SH)
311   CONTINUE

      RETURN
      END
