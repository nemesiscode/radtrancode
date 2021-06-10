      real function lineav(rad,xk)
C     *************************************************************
C     Subroutine to perform a central meridian line-averaged of a Disc
C     using the Minnaert-approximation. Code needs a precomputed file that
C     lists the line-average I/F as a function of minnaert-k, 
C     assuming I(zen=0)=0.0
C
C     Input parameters
C	rad	real	Radiance of I/F at zenith angle=0.0
C	xk	real	Minnaert k-value
C
C     Output parameters
C	lineav	real	line-integrated I/F or radiance
C
C     Pat Irwin		Original 	May 2021
C			Documented	10/6/21
C
C     *************************************************************

      real rad,xk
      integer mtab,ntab
      parameter (mtab=121)
      real Itab(mtab),ktab(mtab)
      common /lineminnaert/ntab,Itab,ktab
 
      if(ntab.ne.121)then
       print*,ntab,' Reading lineminnaert.txt'
       open(12,file='lineminnaert.txt')
        read(12,*)ntab
        do i=1,ntab
         read(12,*)ktab(i),Itab(i)
        enddo
       close(12)
      endif

      i1 = 1+int(xk*100)
      if(i1.ge.ntab)then
       i1=ntab-1
      endif
      frac = 100*(xk-ktab(i1))

C      print*,i1
C      print*,ktab(i1),Itab(i1)
C      print*,ktab(i1+1),Itab(i1+1)
      x = (1.0-frac)*Itab(i1)+frac*Itab(i1+1)      
C      print*,xk,frac,x

      lineav = rad*x

      return

      end
