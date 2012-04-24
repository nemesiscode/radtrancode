      subroutine eigenerr(nx,st,err)
C     $Id:
C     *****************************************************************
C     Subroutine using Numerical Recipes algorithms to analyse covariance
C     matrix in terms of its eigenvectors in order to gain a better
C     estimate of the retrieved error. 
C
C     See Houghton,Taylor,Rodgers, p 133
C
C     Input variables
C	nx	integer	Number of elements in state vector
C	st(mx,mx) real	Covariance matrix of retrieved state vector
C
C     Output variables
C	err(mx)	real	Retrieved errors
C
C     Pat Irwin	2/5/01		Original
C     Pat Irwin	17/10/03	Tidied for Nemesis
C
C     *****************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer nx,i,j
      real st(mx,mx),s1(mx,mx),D(mx),E(mx),err(mx),sum

C     Transfer st to temporary array s1
      do i=1,nx
       do j=1,nx
        s1(i,j)=st(i,j)
       enddo
      enddo

C     Make sure that s1 really is symmetric
      do i=1,nx
       do j=i,nx
        s1(i,j)=st(j,i)
       enddo
      enddo

C     Householder routine to reduce matrix to tridiagonal form 
      call tred2(s1,nx,mx,D,E)

C     Calculate the eigenvalues (D) and eigenvectors (s1)
      call tqli(D,E,nx,mx,s1)

C     Normalise the eigenfunctions
      do j=1,nx
       sum=0.0
       do i=1,nx
        sum=sum+s1(i,j)**2
       enddo
       do i=1,nx
        s1(i,j)=s1(i,j)/sqrt(sum)
       enddo
      enddo

C     Now calculate error
      do i=1,nx
       err(i)=0.0
       do j=1,nx
        err(i)=err(i)+abs(D(j)*s1(i,j)**2)
       enddo
       err(i)=sqrt(err(i))
      enddo

      return

      end      
 
