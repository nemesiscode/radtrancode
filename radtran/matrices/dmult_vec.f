      SUBROUTINE dmult_vec(idim,a,a1,a2,b,b1,c,c1)
C     $Id: dmult_vec.f,v 1.1 2006/08/03 09:18:39 irwin Exp $
C***********************************************************************
C_TITL:	DMULT_VEC.f
C
C_DESC:	This subroutine performs the matrix multiplication C=A*B
C	where	the matrix A has a1 rows and a2 columns
C		the vector B has b1 rows
C		the vector C has c1 rows
C
C_ARGS:	Input variables:
C	idim		INTEGER			Maximum number of elements
C						in the input matrices.
C	a(idim,idim)	DOUBLE PRECISION	Input matrix A.
C	a1		INTEGER			Number of rows in A.
C	a2		INTEGER			Number of columns in A.
C	b(idim,1)	DOUBLE PRECISION	Input matrix B.
C	b1		INTEGER			Number of rows in B.
C
C	Ouput variables:
C	c(idim,1)	DOUBLE PRECISION	Output matrix C.
C	c1		INTEGER			Number of rows in C.
C
C_FILE:	No files openned.
C
C_CALL:	No calls made.
C
C_HIST:	19/1/94	PGJI
C	3jul97	AL	Added cik variable for efficiency.
C***************************** VARIABLES *******************************

      INTEGER i,s
      INTEGER idim,a1,a2,b1,c1
      DOUBLE PRECISION a(idim,idim),b(idim,1),c(idim,1),cik

      c1 = a1
      DO i=1,c1
        cik = 0.0D0
        DO s=1,a2
            cik = cik + a(i,s)*b(s,1)
        ENDDO
        c(i,1) = cik
      ENDDO

      RETURN

      END
C***********************************************************************
C***********************************************************************
