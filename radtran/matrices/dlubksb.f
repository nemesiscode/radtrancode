      SUBROUTINE DLUBKSB(A,NDIM,MDIM,INDX,B)
C***********************************************************************
C_TITL:	DLUBKSB.f
C
C_DESC:	Solves a set of NDIM linear equations Ax=B (in double precision)
C	by both forward and back substitution. This routine takes into
C	account the possibility that B will begin with many zero elements, 
C	so it is efficient for use in matrix inversion.
C
C	The routine is a Numerical Recipies subroutine and is copied
C	almost exactly from the book.
C
C_ARGS: Input variables:
C	A(MDIM,MDIM)	DOUBLE PRECISION	LU decomposition from
C						DLUDCMP.f.
C	NDIM		INTEGER			Number of linear
C						equations.
C	MDIM		INTEGER			Input matrix dimensions.
C	INDX		INTEGER			Input permutation vector
C						from DLUDCMP.f.
C
C       Output variable:
C	B(NDIM)		DOUBLE PRECISION	Solution vector x.
C
C_FILE:	No files openned.
C
C_CALL:	No calls made.
C
C_HIST: ORIGINAL VERSION
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INTEGER I,J,II,LL,NDIM,MDIM
      INTEGER INDX(MDIM)
      DOUBLE PRECISION A(MDIM,MDIM),B(MDIM),SUM

C******************************** CODE *********************************

C When II is set to a positive value, it will become the index of the
C first nonvanishing element of B.
      II = 0

C Do the forward substitution.
      DO 12 I=1,NDIM
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0)THEN
          DO 11 J=II,I-1
            SUM = SUM - A(I,J)*B(J)
11        CONTINUE
        ELSE IF(SUM.NE.0.)THEN
          II = I
        ENDIF
        B(I) = SUM
12    CONTINUE

C Now do the backsubstitution.
      DO 14 I=NDIM,1,-1
        SUM = B(I)
        IF(I.LT.NDIM)THEN
          DO 13 J=I+1,NDIM
            SUM = SUM - A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I) = SUM/A(I,I)
14    CONTINUE

      RETURN

      END
C***********************************************************************
C***********************************************************************
