      real function sproduct(a,b)
C     *******************************************************
C     Perform scalar product of two 3-D vectors a and b.
C     Input variables
C	a(3)	real	vector a
C	b(3)	real	vector b
C
C     Output variable: a.b
C
C     Pat Irwin	Original	20/3/15
C
C     *******************************************************
      real a(3),b(3),sum
      integer i

      sum=0.
      do i=1,3 
       sum=sum+a(i)*b(i)
      enddo
 
      sproduct=sum

      return

      end
