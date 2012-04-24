      SUBROUTINE DINVMARQ(JMOD,ICHECK,A,NDIM,MDIM,B,ERR)
C       $Id:
C***********************************************************************
C_TITL:	MATINV.f
C
C_DESC:	Matrix inversion routine. On exit, AINV contains the MDIMxMDIM
C	inverse of A, and B contains the linear solution
C
C	The routine uses Numerical Recipies subroutines and is copied
C	almost exactly from the book.
C
C_ARGS:	Input variables:
C	A(MDIM,MDIM)	DOUBLE		Input matrix.
C	B(MDIM,1)	DOUBLE		Input vector
C	NDIM		INTEGER		Output matrix dimensions.
C	MDIM		INTEGER		Input matrix dimensions.
C	JMOD		INTEGER		0 = LU, 1 = Gaussian method
C	ICHECK		INTEGER		Set to non-zero to check result
C
C	Output variable:
C	A(MDIM,MDIM)	DOUBLE		Output matrix/Input matrix
C					inversed.
C	B(MDIM,1)	DOUBLE		Linear solution X
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
      PARAMETER (TINY=1E-6)
      INTEGER I,J,NDIM,MDIM,N1,N2,M,ICHECK,JMOD
      INTEGER INDX(IDIM),ERR,I1,J1
C INDX: Index of the row permutation done within LUDCMP.f.
      DOUBLE PRECISION A(MDIM,MDIM),AINV(IDIM,IDIM),IDE(IDIM,IDIM)
      DOUBLE PRECISION D,B(MDIM,1),A1(IDIM,IDIM),DX,B1(IDIM,1)

C******************************** CODE *********************************

      IF(NDIM.GT.MDIM)THEN
        WRITE(*,*)'DINVMARQ :: *ERROR* NDIM > MDIM. Stopping  program.'
        WRITE(*,*)' '
        WRITE(*,*)'NDIM, MDIM = ',ndim,mdim
        STOP
      ENDIF

      IF(MDIM.GT.IDIM)THEN
        WRITE(*,*)'DINVMARQ :: *ERROR* MDIM > IDIM. Stopping  program.'
        WRITE(*,*)' '
        WRITE(*,*)'MDIM, IDIM = ',mdim,idim
        STOP
      ENDIF


      IF(JMOD.EQ.0)THEN
C         LU - DECOMPOSITION
C          print*,'LU decomposition'
C         transfer input array and set up an identity matrix.
          DO 20 J=1,NDIM
           DO 30 I=1,NDIM
            A1(I,J) = A(I,J)
            AINV(I,J) = 0.0
30         CONTINUE
           AINV(J,J) = 1.0
           B1(J,1) = B(J,1)
20        CONTINUE

C         Decompose the matrix just once.
          CALL DLUDCMP1(A1,NDIM,MDIM,INDX,D,ERR)
          IF(ERR.EQ.-1)THEN
           PRINT*,'Error on dinvmarq: aborting'
           RETURN
          ENDIF
C         Find the inverse by columns: it is necessary to recognise that 
C         FORTRAN stores two dimensional data matrices by column, so that 
C         AINV(1,J) is the address of the jth column of AINV.
          DO 40 J=1,NDIM
           CALL DLUBKSB(A1,NDIM,MDIM,INDX,AINV(1,J))
40        CONTINUE

          CALL DMULT_VEC(MDIM,AINV,NDIM,NDIM,B1,NDIM,B,NDIM)

      ELSE
          
C          print*,'Gaussian elimination'
C         transfer input array and set up an identity matrix.
          DO 21 J=1,NDIM
           DO 31 I=1,NDIM
            A1(I,J) = A(I,J)
31         CONTINUE
21        CONTINUE

          M=1
          CALL GAUSSJD(A1,NDIM,MDIM,B,M,M,ERR)

          IF (ERR.EQ.-1)THEN
           PRINT*,'Error in DINVMARQ - GAUSS METHOD HAS COLLAPSED'
           RETURN
          ENDIF

          DO 51 I = 1,NDIM
           DO 61 J = 1,NDIM
            AINV(I,J)=A1(I,J)
61         CONTINUE
51        CONTINUE         

      ENDIF

      IF(ICHECK.NE.0)THEN

       CALL DMULT_MAT(MDIM,A,NDIM,NDIM,AINV,NDIM,NDIM,IDE,N1,N2)
       
       DO 80 I=1,NDIM
        DO 85 J=1,NDIM
         IF(I.EQ.J)THEN
          DX = ABS(1.0-IDE(I,J))
         ELSE
          DX = ABS(IDE(I,J))
         ENDIF

         ERR = 0
         IF(DX.GT.TINY)THEN
C          PRINT*,'ERROR IN DINVERTM'
C          PRINT*,'INVERSE MATRIX DOESNT WORK'
C          PRINT*,'NDIM ',NDIM
C          PRINT*,'DX, TINY = ',DX,TINY

          ERR = -1

C          OPEN(34,FILE='matrices.dat',STATUS='unknown')
C          DO I1=1,NDIM
C           WRITE(34,*)(A(I1,J1),J1=1,NDIM)
C          ENDDO
C          DO I1=1,NDIM
C           WRITE(34,*)(AINV(I1,J1),J1=1,NDIM)
C          ENDDO
C          DO I1=1,NDIM
C           WRITE(34,*)(IDE(I1,J1),J1=1,NDIM)
C          ENDDO
C          CLOSE(34)
          GOTO 101
         ENDIF

85      CONTINUE
80     CONTINUE

C      PRINT*,'Matrix inversion OK'
101    CONTINUE
      ENDIF

      DO 71 I=1,NDIM
       DO 72 J=1,NDIM
        A(J,I)=AINV(J,I)
72     CONTINUE
71    CONTINUE

      RETURN

      END
C***********************************************************************
C***********************************************************************
