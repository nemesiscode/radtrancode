      SUBROUTINE dmult_mat(idim,a,a1,a2,b,b1,b2,c,c1,c2)
C     $Id: dmult_mat.f,v 1.2 2003/04/28 17:08:32 parrish Exp $
C***********************************************************************
C_TITL:	DMULT_MAT.f
C
C_DESC:	This subroutine performs the matrix multiplication C=A*B
C	where	the matrix A has a1 rows and a2 columns
C		the matrix B has b1 rows and b2 columns
C		the matrix C has c1 rows and c2 columns
C
C_ARGS:	Input variables:
C	idim		INTEGER			Maximum number of elements
C						in the input matrices.
C	a(idim,idim)	DOUBLE PRECISION	Input matrix A.
C	a1		INTEGER			Number of rows in A.
C	a2		INTEGER			Number of columns in A.
C	b(idim,idim)	DOUBLE PRECISION	Input matrix B.
C	b1		INTEGER			Number of rows in B.
C	b2		INTEGER			Number of columns in B.
C
C	Ouput variables:
C	c(idim,idim)	DOUBLE PRECISION	Output matrix C.
C	c1		INTEGER			Number of rows in C.
C	c2		INTEGER			Number of columns in C.
C
C_FILE:	No files openned.
C
C_CALL:	No calls made.
C
C_HIST:	19/1/94	PGJI
C	3jul97	AL	Added cik variable for efficiency.
C***************************** VARIABLES *******************************

      INTEGER i,k,s
      INTEGER idim,a1,a2,b1,b2,c1,c2
      DOUBLE PRECISION a(idim,idim),b(idim,idim),c(idim,idim),cik

C First check that the colums of A = rows of B
cc      IF(a2.NE.b1)THEN
cc        WRITE(*,*)'DMULT_MAT.f :: *ERROR* Columns A <> Rows B'
cc        WRITE(*,*)'Stopping program.'
cc        WRITE(*,*)' '
cc        WRITE(*,*)'colums of A, rows of B = ',a2,b1
cc        STOP
cc      ENDIF

      c1 = a1
      c2 = b2
      DO k=1,c2
         DO i=1,c1
            cik = 0.0D0
            DO s=1,a2
               cik = cik + a(i,s)*b(s,k)
            ENDDO
            c(i,k) = cik
         ENDDO
      ENDDO

      RETURN

      END
C***********************************************************************
C***********************************************************************
