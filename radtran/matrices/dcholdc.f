      SUBROUTINE DCHOLDC(A,NDIM,MDIM,P)
C***********************************************************************
C_TITL: DLUDCMP.f
C
C_DESC: Given a NDIMxNDIM matrix A, with physical dimension MDIM, this
C       routine constructs its Cholesky decomposition, A=L.L^T. 
C       On input, only the upper triangle of A is needed. It is not
C       modified. The Cholesky factor L is returned in the lower triangle
C       of A, except for its diagonal components which are returned in
C       P(1:NDIM)
C
C       The routine is a Numerical Recipies subroutine and is copied 
C       almost exactly from the book.
C
C_ARGS: Input variables:
C       A(MDIM,MDIM)    DOUBLE PRECISION        Input matrix. Only
C						Upper triangle form is
C						needed. Lower triangle is
C						overwritten.	
C       NDIM            INTEGER                 Number of elements in
C                                               input matrix.
C       MDIM            INTEGER                 Maximum number of 
C                                               elements in input matrix.
C
C       Output variables:
C       P(NDIM)         DOUBLE PRECISION        Diagonal elements of
C						Cholesky factor L. 
C	A(NDIM,NDIM) 	DOUBLE PRECISION	Lower traingle gets 
C						overwritten by Cholesky
C						factor L.
C
C_FILE: No files opened.
C
C_CALL: No calls made.
C
C_HIST: 22nov04 PGJI     ORIGINAL VERSION
C***************************** VARIABLES *******************************
      IMPLICIT NONE
C     Read in IDIM from common file.
      INCLUDE '../includes/arrdef.f'
      INTEGER NDIM,MDIM
      DOUBLE PRECISION A(MDIM,MDIM),P(MDIM),SUM,A1(IDIM,IDIM)
      INTEGER i,j,k,l
 
      IF(NDIM.GT.MDIM)THEN
       PRINT*,'Error: DCHOLDC.F: NDIM > MDIM',NDIM,MDIM
       STOP
      ENDIF

      IF(MDIM.GT.IDIM)THEN
       PRINT*,'Error: DCHOLDC.F: MDIM > IDIM',MDIM,IDIM
       STOP
      ENDIF

C     Transfer to read-only array in case of abort.
      do i=1,ndim
       do j=1,ndim
        a1(i,j)=a(i,j)
       enddo
      enddo

      do 13 i=1,ndim
        do 12 j=i,ndim
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(sum.le.0.)then
             print*,'dcholdc failed. Matrix A is not positive-definite'
             open(12,file='dcholdc.dat',status='unknown')
              write(12,*)ndim
              do k=1,ndim
               write(12,*)(a1(i,k),l=1,ndim)
              enddo
             close(12)
             stop
            endif		
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
12      continue
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.
