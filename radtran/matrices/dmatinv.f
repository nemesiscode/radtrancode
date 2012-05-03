      SUBROUTINE DMATINV(A,NDIM,MDIM,AINV)
C***********************************************************************
C_TITL:	DMATINV.f
C
C_DESC:	Matrix inversion routine in double precision. On exit, AINV 
C	contains the MDIMxMDIM inverse of A where the original A has
C	been destroyed.
C
C	The routine used Numerical Recipies subroutines and is copied
C	almost exactly from the book.
C
C_ARGS:	Input variables:
C	A(MDIM,MDIM)	DOUBLE PRECISION	Input matrix.
C	NDIM		INTEGER			Number of elements in
C						input matrix.
C	MDIM		INTEGER			Maximum number of
C						elements in input matrix.
C
C	Output variable:
C	AINV(MDIM,MDIM)	DOUBLE PRECISION	Matrix A inversed.
C
C_FILE:	No files openned.
C
C_CALL:	dludcmp	Replaces matrix A with a LU (lower/upper triangular)
C		decomposition of a rowwise permutation of itself (in
C		double precision).
C	dlubksb	Solves a set of NDIM linear equations Ax=B by both 
C		forward and back substitution (in double precision).
C
C_HIST:	21jan88	SBC	ORIGINAL VERSION
C	15apr93	PGJI	Modified
C***************************** VARIABLES *******************************

      IMPLICIT NONE
C     Read in IDIM from common file.
      INCLUDE '../includes/arrdef.f'
      INTEGER I,J,NDIM,MDIM
      INTEGER INDX(IDIM)
C INDX: Index of the row permutation done within DLUDCMP.f.
      DOUBLE PRECISION D,A(MDIM,MDIM),AINV(MDIM,MDIM)

C******************************** CODE *********************************

      IF(NDIM.GT.MDIM)THEN
        WRITE(*,*)'DMATINV.f :: *ERROR* NDIM > MDIM. Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'NDIM, MDIM = ',ndim,mdim
        STOP
      ENDIF
      IF(MDIM.GT.IDIM)THEN
        WRITE(*,*)'DMATINV.f :: *ERROR* MDIM > IDIM. Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'MDIM, IDIM = ',Mdim,idim
        STOP
      ENDIF
C Set up an identity matrix.
      DO 20 J=1,NDIM
        DO 30 I=1,NDIM
          AINV(I,J) = 0.0
30      CONTINUE
        AINV(J,J) = 1.0
20    CONTINUE

C LU decompose the matrix (just the once).
      CALL DLUDCMP(A,NDIM,MDIM,INDX,D)

C Find the inverse by columns: it is necessary to recognise that FORTRAN  
C stores two dimensional data matrices by column, so that AINV(1,J) is the
C address of the jth column of AINV.
      DO 40 J=1,NDIM
        CALL DLUBKSB(A,NDIM,MDIM,INDX,AINV(1,J))
40    CONTINUE

      RETURN

      END
C***********************************************************************
C***********************************************************************
