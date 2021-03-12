      SUBROUTINE DINVERTM(JMOD,ICHECK,A,NDIM,MDIM,AINV)
C       $Id:
C***********************************************************************
C_TITL:	MATINV.f
C
C_DESC:	Matrix inversion routine. On exit, AINV contains the MDIMxMDIM
C	inverse of A.
C
C	The routine uses Numerical Recipies subroutines and is copied
C	almost exactly from the book.
C
C_ARGS:	Input variables:
C	A(MDIM,MDIM)	DOUBLE		Input matrix.
C	NDIM		INTEGER		Output matrix dimensions.
C	MDIM		INTEGER		Input matrix dimensions.
C	JMOD		INTEGER		0 = LU, 1 = Gaussian method
C					2 = Cholesky
C
C	Output variable:
C	AINV(MDIM,MDIM)	DOUBLE		Output matrix/Input matrix
C					inversed.
C	ICHECK		INTEGER		0 if inversion OK, 1 if not.
C
C_FILE:	No files openned.
C
C_CALL:	ludcmp	Replaces matrix A with a LU (lower/upper triangular)
C		decomposition of a rowwise permutation of itself.
C	lubksb	Solves a set of NDIM linear equations Ax=B by both
C		forward and back substitution.
C
C_HIST:	21jan88 SBC	ORIGINAL VERSION
C	15apr93 PGJI	Modified
C***************************** VARIABLES *******************************

      IMPLICIT NONE
C     Read in IDIM from common file.
      INCLUDE '../includes/arrdef.f'
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1E-3)
      INTEGER I,J,NDIM,MDIM,N1,N2,M,ICHECK,JMOD,I1,J1
      INTEGER INDX(IDIM)
C     INDX: Index of the row permutation done within LUDCMP.f.
      DOUBLE PRECISION A(MDIM,MDIM),AINV(MDIM,MDIM),IDE(IDIM,IDIM)
      DOUBLE PRECISION D,B(IDIM,1),A1(IDIM,IDIM),DX,P(IDIM),X(IDIM)
      DOUBLE PRECISION B1(IDIM),DXMAX,AINV1(IDIM,IDIM)
      LOGICAL NTEST,ISNAN
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C******************************** CODE *********************************

      ICHECK=0
      IF(NDIM.GT.MDIM)THEN
        WRITE(*,*)'DINVERTM :: *ERROR* NDIM > MDIM. Stopping  program.'
        WRITE(*,*)' '
        WRITE(*,*)'NDIM, MDIM = ',ndim,mdim
        STOP
      ENDIF

      IF(MDIM.GT.IDIM)THEN
        WRITE(*,*)'DINVERTM :: *ERROR* MDIM > IDIM. Stopping  program.'
        WRITE(*,*)' '
        WRITE(*,*)'MDIM, IDIM = ',mdim,idim
        STOP
      ENDIF

      IF(JMOD.EQ.0)THEN
          if(idiag.gt.0)print*,'dinvertm : LU decomposition'
C         LU - DECOMPOSITION
C         transfer input array and set up an identity matrix.
          DO 20 J=1,NDIM
           DO 30 I=1,NDIM
            A1(I,J) = A(I,J)
            AINV(I,J) = 0.0
30         CONTINUE
           AINV(J,J) = 1.0
20        CONTINUE

C         Decompose the matrix just once.
          CALL DLUDCMP(A1,NDIM,IDIM,INDX,D)
C         Find the inverse by columns: it is necessary to recognise that 
C         FORTRAN stores two dimensional data matrices by column, so that 
C         AINV(1,J) is the address of the jth column of AINV.
          DO 40 J=1,NDIM
           CALL DLUBKSB(A1,NDIM,IDIM,INDX,AINV1(1,J))
40        CONTINUE

          DO I=1,NDIM
           DO J=1,NDIM
            AINV(I,J)=AINV1(I,J)
           ENDDO
          ENDDO

      ELSEIF(JMOD.EQ.1)THEN
          
          if(idiag.gt.0)print*,'dinvertm: Gaussian elimination'
C         transfer input array and set up an identity matrix.
          DO 21 J=1,NDIM
           DO 31 I=1,NDIM
            A1(I,J) = A(I,J)
31         CONTINUE
           B(J,1) = 1.0
21        CONTINUE

          M=1
          CALL DGAUSSJ(A1,NDIM,IDIM,B,M,M)

          DO 51 I = 1,NDIM
           DO 61 J = 1,NDIM
            AINV(I,J)=A1(I,J)
61         CONTINUE
51        CONTINUE         

      ELSE
C         CHOLESKY DECOMPOSITION
          if(idiag.gt.0)print*,'dinvertm: Cholesky decomposition'
C         transfer input array and to upper triangle of A1 and set up 
C         an identity matrix.
          DO 25 J=1,NDIM
           DO 34 I=1,J-1
            AINV(I,J)=0.0
34         CONTINUE
           DO 35 I=J,NDIM
            A1(J,I) = A(J,I)
            AINV(I,J) = 0.0
35         CONTINUE
           AINV(J,J) = 1.0
25        CONTINUE

C         Decompose the matrix just once.
          CALL DCHOLDC(A1,NDIM,IDIM,P)

          DO 45 J=1,NDIM
           DO 44 I=1,NDIM
            B1(I)=AINV(I,J)
44         CONTINUE
           CALL DCHOLSL(A1,NDIM,IDIM,P,B1,X)
           DO 46 I=1,NDIM
            AINV(I,J)=X(I)
46         CONTINUE
45        CONTINUE

      ENDIF

C     Checking inversion is OK

      DO I=1,NDIM
       DO J=1,NDIM
        A1(I,J)=A(I,J)
        AINV1(I,J)=AINV(I,J)
        NTEST=ISNAN(AINV1(I,J))
        IF(NTEST)ICHECK=1
       ENDDO
      ENDDO

      CALL DMULT_MAT(IDIM,A1,NDIM,NDIM,AINV1,NDIM,NDIM,IDE,N1,N2)
       
      DO 80 I=1,NDIM
        DO 85 J=1,NDIM
         IF(I.EQ.J)THEN
          DX = ABS(1.0-IDE(I,J))
         ELSE
          DX = ABS(IDE(I,J))
         ENDIF
         IF(DX.GT.TINY)THEN
          ICHECK = 1
          DXMAX = DX
         ENDIF
85      CONTINUE
80    CONTINUE

      IF(ICHECK.EQ.0) THEN
        if(idiag.gt.0)print*,'dinvertm.f: Matrix inversion OK'
      ELSE
        if(idiag.gt.0)print*,'dinvertm.f: Matrix not inverted well'
        if(idiag.gt.0)print*,'DXMAX,TINY',DXMAX,TINY
        if(idiag.gt.0)print*,'JMOD : ',JMOD
      ENDIF

      RETURN

      END
C***********************************************************************
C***********************************************************************
