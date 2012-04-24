      SUBROUTINE DCHOLSL(A,NDIM,MDIM,P,B,X)
C***********************************************************************
C_TITL: DCHOLSL.f
C
C_DESC: Solves a set of NDIM linear equations Ax=B (in double precision)
C       where A is a positive definite symmetric matrix with physical
C	dimension MDIM. A and P are input as the output of the routine
C	DCHOLDC. Only the loer triangle of A is accessed. B(1:NDIM) is 
C	input as the right-hand side vector. The solution vector is
C	returned in x(1:n). A, NDIM, MDIM and P are not modified and
C       can be left in place for successive calls with different
C	right-hand sides B. B is not modified unless you identify B and X
C	in the calling sequence, which is allowed. 
   
C       The routine is a Numerical Recipies subroutine and is copied
C       almost exactly from the book.
C
C_ARGS: Input variables:
C       A(MDIM,MDIM)    DOUBLE PRECISION        Upper triangle form 
C						contains real symmetric 
C						matrix. Lower triangle
C						contains Cholesky factor L.
C       NDIM            INTEGER                 Number of elements in
C                                               input matrix.
C       MDIM            INTEGER                 Maximum number of 
C       P(NDIM)         DOUBLE PRECISION        Diagonal elements of
C                                               Cholesky factor L. 
C
C	B(NDIM)		DOUBLE PRECISION	Right-hand side vector
C
C       Output variable:
C	X(NDIM)		DOUBLE PRECISION	Solution.
C
C_FILE: No files opened.
C                                               
C_CALL: No calls made.
C       
C_HIST: 22nov04 PGJI     ORIGINAL VERSION
C***********************************************************************
      IMPLICIT NONE
      INTEGER NDIM,MDIM
      DOUBLE PRECISION A(MDIM,MDIM),B(MDIM),P(MDIM),X(MDIM)
      INTEGER i,k
      DOUBLE PRECISION sum

      do 12 i=1,ndim
        sum=b(i)
        do 11 k=i-1,1,-1
          sum=sum-a(i,k)*x(k)
11      continue
        x(i)=sum/p(i)
12    continue
      do 14 i=ndim,1,-1
        sum=x(i)
        do 13 k=i+1,ndim
          sum=sum-a(k,i)*x(k)
13      continue
        x(i)=sum/p(i)
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.
