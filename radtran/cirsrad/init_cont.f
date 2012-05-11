      SUBROUTINE INIT_CONT(VMIN,VMAX,WING)
C     ************************************************************
C     Initialise continuum
C
C     ************************************************************
      IMPLICIT NONE

      INTEGER I,J,K,LAYER,L
      REAL VMIN,VMAX,WING
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/contdef.f'
      DOUBLE PRECISION DMATRIX(IORDP1,IORDP1),DUNIT(IORDP1,IORDP1)

      NBIN=1+INT((VMAX-VMIN)/WING)

      IF(NBIN.GT.MAXBIN)THEN       
        WRITE(*,*)'INIT_CONT.f :: *ERROR* NBIN > MAXBIN'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'NBIN, MAXBIN = ',NBIN,MAXBIN
        STOP
      ENDIF

      DO 127 I=1,NBIN
        VBIN(I) = VMIN + (I - 1)*WING
127   CONTINUE

C     Initialising continuum polynomial
      DO 16 I=1,IORDP1
        DO 14 J=1,NBIN
         DO 11 LAYER=1,MAXLAY
           CONTINK(I,LAYER,J) = 0.0       
           CONVALS(I,LAYER,J) = 0.0       
11       CONTINUE
14      CONTINUE
16    CONTINUE

      DO 444 K=1,IORDP1
        CONWAV(K) = FLOAT(K - 1)*WING/FLOAT(IORDER)
444   CONTINUE      

C     Setting up matrix of wavenumbers for use in polynomial calculation
      DO 18 K=1,IORDP1
        MATRIX(1,K) = 1.0
        DMATRIX(1,K) = 1.0
        DO 19 J=2,IORDP1
          MATRIX(J,K) = MATRIX(J-1,K)*CONWAV(K)
          DMATRIX(J,K) = DMATRIX(J-1,K)*CONWAV(K)
19      CONTINUE
18    CONTINUE

      L = IORDP1       
C     Find the inverse of the matrix
      CALL DMATINV(DMATRIX,L,IORDP1,DUNIT)

      DO K=1,IORDP1
       DO J=1,IORDP1
        UNIT(J,K)=SNGL(DUNIT(J,K))
       ENDDO
      ENDDO

      RETURN

      END

