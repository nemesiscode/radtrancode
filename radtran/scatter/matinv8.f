	SUBROUTINE MATINV8(A,N,NDIM,AINV)
C-----------------------------------------------------------------------
C_TITLE: MATINV8: matrix inversion routine
C
C_ARGS:  A id NxN matrix stored in a NDIMxNDIM array. AINV id also a NDIMxNDIM
C        array.
C
C_KEYS:
C
C_DESCR: on exit AINV contains the NxN inverse of A. A is destroyed.
C        The routine used Numerical Recipies subroutines and is copied almost
C        exactly from the book.
C
C_FILES: none
C
C_CALLS: none
C
C_BUGS:
C
C_HIST:  21jan88 SBC ORIGINAL VERSION
C         3may88 LWK: converted to DOUBLE PRECISION
C        30may88 LWK: renamed MATINV8 to avoid conflict with NAG$SHARE module
C
C_END:
C-----------------------------------------------------------------------
	IMPLICIT NONE
        INCLUDE '../includes/arrdef.f'
C	uses an internally declared array to save the bother of remembering
C	how to use Numerical Recipies routines
	INTEGER N,NDIM,INDX(MAXMU),I,J
	DOUBLE PRECISION A(NDIM,NDIM),AINV(NDIM,NDIM),D
	IF(NDIM.GT.MAXMU)THEN
	  WRITE(6,10)
10	  FORMAT(' ERROR - internal array in MATINV8 is too small')
	  STOP
	  END IF
	DO 20 J=1,N
	DO 30 I=1,N
	AINV(I,J)=0.D0
30	CONTINUE
	AINV(J,J)=1.D0
20	CONTINUE
C
	CALL LUDCMP8(A,N,NDIM,INDX,D)
	DO 40 J=1,N
	CALL LUBKSB8(A,N,NDIM,INDX,AINV(1,J))
40	CONTINUE
	RETURN
	END
      SUBROUTINE LUDCMP8(A,N,NP,INDX,D)
C  3/5/88 ...LWK... FROM [ATMRJW.RECIPES], CONVERTED TO DOUBLE PRECISION
      IMPLICIT NONE
      INTEGER NMAX,I,J,K,N,NP,INDX(N),IMAX
      PARAMETER (NMAX=100)
      DOUBLE PRECISION A(NP,NP),VV(NMAX),D,AAMAX,SUM,DUM,TINY
      PARAMETER (TINY=1.0D-20)
      IF(N.GT.NMAX)THEN
       PRINT*,'Error in matinv8:ludcmp8. N>NMAX'
       PRINT*,N,NMAX
       STOP 
      ENDIF
      D=1.0D0
      DO 12 I=1,N
        AAMAX=0.D0
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) THEN
          PRINT*, 'Singular matrix.'
	  STOP
        ENDIF
        VV(I)=1.D0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.D0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1.D0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.D0)A(N,N)=TINY
      RETURN
      END
      SUBROUTINE LUBKSB8(A,N,NP,INDX,B)
C  3/5/88 ...LWK...  FROM [ATMRJW.RECIPES], CONVERTED TO DOUBLE PRECISION
      IMPLICIT NONE
      INTEGER N,NP,INDX(N),I,J,II,LL
      DOUBLE PRECISION A(NP,NP),B(N),SUM
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
