      SUBROUTINE DLUDCMP1(A,NDIM,MDIM,INDX,D,IERR)
C***********************************************************************
C_TITL:	DLUDCMP.f
C
C_DESC:	Given a NDIMxNDIM matrix A, with physical dimension MDIM, this
C	routine replaces it by the LU (lower/upper triangular)
C	decomposition of a rowwise permutation of itself (in double
C	precision). This routine is used in combination with DLUBKSB to
C	solve linear equations or invert a matrix.
C
C	The routine is a Numerical Recipies subroutine and is copied 
C	almost exactly from the book.
C
C_ARGS:	Input variables:
C	A(MDIM,MDIM)	DOUBLE PRECISION	Input matrix.
C	NDIM		INTEGER			Number of elements in
C						input matrix.
C	MDIM		INTEGER			Maximum number of 
C                                               elements in input matrix.
C
C	Output variables:
C	INDX(MDIM)	INTEGER			Record of the row
C						permutation effected by
C						the partial pivoting.
C	D		DOUBLE PRECISION	Output equal to +/-1
C						depending on whether the
C						number of row interchanges
C						was even or odd,
C						respectively.
C	IERR		INTEGER			0 = OK, 1= Error
C
C_FILE: No files openned.
C
C_CALL: No calls made.
C
C_HIST:	ORIGINAL VERSION
C***************************** VARIABLES *******************************

      IMPLICIT NONE
C     Read in IDIM from common file.
      INCLUDE '../includes/arrdef.f'
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1.0E-30)
C TINY: A small number.

      INTEGER I,J,K,IMAX,NDIM,MDIM,INDX(MDIM),IERR
      DOUBLE PRECISION A(MDIM,MDIM),VV(IDIM)
C VV: Stores the implicit scaling of each row.
      DOUBLE PRECISION D,AAMAX,SUM,DUM
C DUM: Figure of merit for the pivot.

C******************************** CODE *********************************

      IF(NDIM.GT.MDIM)THEN
       PRINT*,'Error: DLUDCMP1.F: NDIM > MDIM',NDIM,MDIM
       STOP
      ENDIF

      IF(MDIM.GT.IDIM)THEN
       PRINT*,'Error: DLUDCMP1.F: MDIM > IDIM',MDIM,IDIM
       STOP
      ENDIF

      D = 1.0                          ! no row interchanges yet.
      IERR = 0
C Loop over rows to get the implicit scaling information.
      DO 12 I=1,NDIM
        AAMAX = 0.0
        DO 11 J=1,NDIM
          IF(ABS(A(I,J)).GT.AAMAX)AAMAX = ABS(A(I,J))
11      CONTINUE
        IF(AAMAX.EQ.0.)THEN
          WRITE(*,*)'DLUDCMP1.f :: *ERROR* Singular matrix. Stopping'
          WRITE(*,*)'program.'
          WRITE(*,*)' '
          WRITE(*,*)'AAMAX = 0.0 ==> the largest element cannot be'
          WRITE(*,*)'non-zero.'
          IERR = -1
          RETURN
        ENDIF
        VV(I) = 1.0/AAMAX              ! save the scaling.
12    CONTINUE

C Loop over columns using Crout's method of decompostion.
      DO 19 J=1,NDIM
        IF(J.GT.1)THEN
          DO 14 I=1,J-1
            SUM = A(I,J)
            IF(I.GT.1)THEN
              DO 13 K=1,I-1
                SUM = SUM - A(I,K)*A(K,J)
13            CONTINUE
              A(I,J) = SUM
            ENDIF
14        CONTINUE
        ENDIF

C Initialise for the search of the largest pivot element.
        AAMAX = 0.0
        DO 16 I=J,NDIM
          SUM = A(I,J)
          IF(J.GT.1)THEN
            DO 15 K=1,J-1
              SUM = SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J) = SUM
          ENDIF
          DUM = VV(I)*ABS(SUM)
C Check to see if the figure of merit is better than the best so far ...
          IF(DUM.GE.AAMAX)THEN
            IMAX = I
            AAMAX = DUM
          ENDIF
16      CONTINUE

C Check to see if the interchange of rows is required.
        IF(J.NE.IMAX)THEN
          DO 17 K=1,NDIM
            DUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = DUM
17        CONTINUE
C Change the parity of D ...
          D = -D
C ... and interchange the scale factor if so.
          VV(IMAX) = VV(J)
        ENDIF
        INDX(J) = IMAX

C Now, divide by the pivot element. If the pivot element is zero then the
C matrix is singular (at least to the precision of the algorithm). For 
C some applications on singular matrices, it is desirable to substitute
C TINY for zero.
        IF(J.NE.NDIM)THEN
          IF(A(J,J).EQ.0.)A(J,J) = TINY
          DUM = 1./A(J,J)
          DO 18 I=J+1,NDIM
            A(I,J) = A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE

      IF(A(NDIM,NDIM).EQ.0.)A(NDIM,NDIM) = TINY

      RETURN

      END
C***********************************************************************
C***********************************************************************
