      SUBROUTINE HYDROSTATP(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,HTAN,
     1 PTAN,SCALE)
C     ****************************************************************
C     Subroutine to rescale the pressures of a H/P/T profile according to
C     the hydrostatic equation above and below a specified altitude
C     given the pressure at that altitude.
C
C     Pat Irwin		27/8/06
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      INTEGER NPRO,IPLANET,JTAN,I
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),LATITUDE,TEMP,SH
      REAL RADIUS,G,MOLWT,R,SCALE(MAXPRO),HTAN,PTAN,DELH,P1(MAXPRO)
      PARAMETER (R=8.3143)
      CHARACTER*8 PNAME


C     **************** Code **************

C     First find level immediately below the reference altitude
      DO I=1,NPRO
       IF(H(I).LT.HTAN)THEN
        JTAN=I
       ENDIF
       P1(I)=P(I)
       CALL NEWGRAV(IPLANET,LATITUDE,H(I),RADIUS,G,PNAME)
       SCALE(I)=R*T(I)/(MOLWT*G)
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

      open(12,file='hydrostatP.dat',status='unknown')
      write(12,*)npro
      write(12,*)'Htan, Ptan = ',HTAN,PTAN
      write(12,*)'H   T   P   New P'
      do i=1,npro
       write(12,*)h(i),t(i),p1(i),p(i)
      enddo
      close(12)

      RETURN
      END
