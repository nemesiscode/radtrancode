      SUBROUTINE HYDROSTATHMOD(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,SCALE)
C     ****************************************************************
C     Subroutine to rescale the heights of a H/P/T profile according to
C     the hydrostatic equation above and below the level where height=0.
C
C     Pat Irwin		27/8/06
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      INTEGER NPRO,IPLANET,JZERO,I,NPRO1
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),LATITUDE,DELH,X,TEMP
      REAL RADIUS,G,MOLWT,R,SCALE(MAXPRO),XMASS(MAXPRO),MOLWT1,PTMP
      PARAMETER (R=8.3143)
      CHARACTER*8 PNAME
      CHARACTER*100 BUFFER

1     FORMAT(A)
      OPEN(12,FILE='molwt.prf',STATUS='old')
      READ(12,*)NPRO1
      IF(NPRO1.NE.NPRO)THEN
       PRINT*,'Error in hydrostatHmod.f NPRO <> NPRO1'
       STOP
      ENDIF
      READ(12,1)BUFFER


C     **************** Code **************
C     First find level closest to zero altitude
      DELH=1000.0
      DO I=1,NPRO
       X = ABS(H(I))
       IF(X.LT.DELH)THEN
        DELH=X
        JZERO=I
       ENDIF
       READ(12,*)PTMP,XMASS(I)
      ENDDO
      CLOSE(12)

      print*,'JZERO = ',JZERO
      IF(JZERO.GT.1.AND.JZERO.LT.NPRO)THEN
       H(JZERO)= 0.0
      ENDIF

      DO 301 I=JZERO+1,NPRO
           TEMP = 0.5*(T(I-1) + T(I))           
           CALL NEWGRAV(IPLANET,LATITUDE,H(I-1),RADIUS,G,PNAME)
           MOLWT1 = 0.5*(XMASS(I-1) + XMASS(I))
           SCALE(I)=R*TEMP/(MOLWT1*G)
           H(I) = H(I-1) - SCALE(I)*LOG(P(I)/P(I-1))
301   CONTINUE

      DO 311 I=JZERO-1,1,-1
           TEMP = 0.5*(T(I) + T(I+1))
           CALL NEWGRAV(IPLANET,LATITUDE,H(I+1),RADIUS,G,PNAME)
           MOLWT1 = 0.5*(XMASS(I) + XMASS(I+1))
           SCALE(I)=R*TEMP/(MOLWT1*G)
           H(I) = H(I+1) - SCALE(I)*LOG(P(I)/P(I+1))
C           print*,I,TEMP,P(I),SCALE(I),H(I)
311   CONTINUE


      IF(JZERO.GT.1.AND.JZERO.LT.NPRO)THEN
        SCALE(JZERO)=0.5*(SCALE(JZERO+1)+SCALE(JZERO-1))
      ENDIF
      IF(JZERO.EQ.NPRO)SCALE(JZERO)=SCALE(JZERO-1)
      IF(JZERO.EQ.1)SCALE(JZERO)=SCALE(JZERO+1)

      RETURN
      END
