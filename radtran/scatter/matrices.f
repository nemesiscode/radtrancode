      SUBROUTINE MEQU(AMAT1,N,IDIM,AMAT2)
C     $Id: matrices.f,v 1.2 2011-06-17 15:57:53 irwin Exp $
C     *************************************************************
C     Subroutine to transfer an NxN matrix contained in AMAT2 into AMAT1
C
C     Input variables
C	AMAT2(IDIM,IDIM)	DOUBLE	Input matrix array
C	N			INTEGER	Size of matrix (NxN)
C	IDIM			INTEGER	Dimension of array
C
C     Output variable:
C	AMAT1(IDIM,IDIM)	DOUBLE	Transfered matrix
C
C     Code commented by Pat Irwin  17/9/96
C
C     *************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AMAT1(IDIM,IDIM),AMAT2(IDIM,IDIM)
      DO J=1,N
	DO I=1,N
	  AMAT1(I,J)=AMAT2(I,J)
	ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE MADD(CONST,AM1,AM2,ANS,N1,N2,IDM1,IDM2)
C     *************************************************************
C     Subroutine to add one matrix to the multiple of another:
C	ANS = AM1 + CONST*AM2
C
C     Input variables:
C	CONST		DOUBLE	Multiplicative constant
C	AM1(IDM1,IDM2)	DOUBLE	1st matrix array
C	AM2(IDM1,IDM2)	DOUBLE	2nd matrix array
C	N1		INTEGER	Number of rows used
C	N2		INTEGER	Number of columns used
C
C      Output variable
C	ANS(IDM1,IDM2)	DOUBLE	Output matrix
C
C     Code commented by Pat Irwin  17/9/96
C
C     *************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AM1(IDM1,IDM2),AM2(IDM1,IDM2),ANS(IDM1,IDM2)
      DO J=1,N2
	DO I=1,N1
	  ANS(I,J)=AM1(I,J)+CONST*AM2(I,J)
	ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE MMUL(CONST,AM1,AM2,ANS,N1,N2,N3,IDM1,IDM2,IDM3)
C     *************************************************************
C     Subroutine to multiply two  matrices together
C	ANS = CONST*AM1*AM2
C
C     Input variables:
C	CONST		DOUBLE	Multiplicative constant
C	AM1(IDM1,IDM2)	DOUBLE	1st matrix array
C	AM2(IDM1,IDM3)	DOUBLE	2nd matrix array
C	N1		INTEGER	Number of rows used of 1st matrix
C	N2		INTEGER	Number of columns used of 1st matrix and number
C				of rows of 2nd matrix
C	N3		INTEGER	Number of columns used of 2nd matrix
C	IDM1		INTEGER No. rows of 1st array
C	IDM2		INTEGER No. cols of 1st array and no. rows of 2nd array
C	IDM3		INTEGER	No. of cols of 2nd array
C
C      Output variable
C	ANS(IDM1,IDM3)	DOUBLE	Output matrix
C
C     Code commented by Pat Irwin  17/9/96
C     Added variable AIJ and moved CONST from K loop for efficiency AL 3-JUL-97
C
C     *************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AM1(IDM1,IDM2),AM2(IDM2,IDM3),ANS(IDM1,IDM3)
      DO J=1,N3
	DO I=1,N1
          AIJ = 0.0D0
	  DO K=1,N2
	    AIJ = AIJ + AM1(I,K)*AM2(K,J)
	  ENDDO
          ANS(I,J)=CONST*AIJ
	ENDDO
      ENDDO
      RETURN
      END

