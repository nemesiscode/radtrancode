      subroutine dadd_mat(idim,a,a1,a2,b,b1,b2,c,c1,c2)
C     $Id: dadd_mat.f,v 1.1.1.1 2000/08/17 09:27:01 irwin Exp $
C     ***********************************************************
C     This subroutine performs the matrix addition C=A+B
C     where the matrix A has a1 rows and a2 columns
C           the matrix B has b1 rows and b2 columns
C           the matrix C has c1 rows and c2 columns
C
C     Pat Irwin   10/9/96
C	ALW	11/7/99
C     ***********************************************************
      integer idim,a1,a2,b1,b2,c1,c2
      double precision a(idim,idim),b(idim,idim),c(idim,idim)
      integer i,k,s

C     first check that the rows of A = rows of B
      if(a1.ne.b1) then
       print*,'Add_mat: Rows A <> Rows B'
       stop
      end if
C     then check that the columns of A = columns of B
      if(a2.ne.b2) then
       print*,'Add_mat: Columns A <> Columns B'
       stop
      end if

C     If OK then continue
      c1=a1
      c2=a2
      do 20 i=1,c1
       do 10 k=1,c2
         c(i,k)=a(i,k)+b(i,k)
10     continue
20    continue
      return
      end       
