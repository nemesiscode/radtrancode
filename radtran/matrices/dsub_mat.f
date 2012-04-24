      subroutine dsub_mat(idim,a,a1,a2,b,b1,b2,c,c1,c2)
C     $Id: dsub_mat.f,v 1.1 2002/02/04 11:28:06 irwin Exp $
C     ***********************************************************
C     This subroutine performs the matrix subtraction C=A-B
C     where the matrix A has a1 rows and a2 columns
C           the matrix B has b1 rows and b2 columns
C           the matrix C has c1 rows and c2 columns
C
C     Pat Irwin   4/10/96
C     ***********************************************************
      integer idim,a1,a2,b1,b2,c1,c2
      double precision a(idim,idim),b(idim,idim),c(idim,idim)
      integer i,k,s
C     first check that the rows of A = rows of B
      if(a1.ne.b1) then
       print*,'Rows A <> Rows B'
       stop
      end if
C     then check that the columns of A = columns of B
      if(a2.ne.b2) then
       print*,'Columns A <> Columns B'
       stop
      end if
      c1=a1
      c2=a2
      do 20 i=1,c1
       do 10 k=1,c2
         c(i,k) = a(i,k) - b(i,k)
10     continue 
20    continue
      return
      end       
