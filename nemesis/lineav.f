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

      real rad,xk,xk1
      integer mtab,ntab
      parameter (mtab=151)
      real Itab(mtab),ktab(mtab)
      logical isnan
      character*100 aname
      common /lineminnaert/ntab,Itab,ktab
 
      if(ntab.ne.151)then
       print*,ntab,'Reading lineminnaert.txt'
       aname='lineminnaert.txt'
       call datarchive(aname) 
       open(12,file=aname,status='old')
        read(12,*)ntab
        print*,'ntab = ',ntab
        do i=1,ntab
         read(12,*)ktab(i),Itab(i)
         print*,i,ktab(i),Itab(i)
        enddo
       close(12)
      endif

C     Add catch to return zero if xk is NAN
      if(isnan(xk))then
       lineav=0.0
       return
      endif

C     Add catch to not extrapolate to k<ktab(1)
      xk1=xk
      if(xk1.lt.ktab(1))xk1=ktab(1)
C      if(xk1.gt.ktab(ntab))xk1=ktab(ntab)
      i1 = 1+int(xk1*100)
      if(i1.ge.ntab)then
       i1=ntab-1
      endif
      frac = 100*(xk1-ktab(i1))
      x = (1.0-frac)*Itab(i1)+frac*Itab(i1+1)      

      lineav = rad*x

      return

      end
