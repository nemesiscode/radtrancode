      subroutine getquad(itype,a,b,n,weight,abscis)
C     *************************************************************
C     Subroutine to return Gaussian quadrature weights and absissas
C
C     Input variables:
C	itype		integer	Quadature identifier
C      				1 - Gauss-Legendre
C				2 - Gauss-Rational
C				3 - Gauss-Laguerre
C				4 - Gauss-Hermite
C	a		double	Minimum of integration range
C	b		double	Maximum of integration range
C	n		integer	Number of ordinates
C
C     Output variables
C	weight(n)	double	Gaussian weights
C	abscis(n)	double	Abscissae
C
C     Pat Irwin		30/9/96
C
C     *************************************************************
      implicit double precision(a-h,o-z)
      double precision a,b,weight(n),abscis(n)
      integer n,ifail,itype,itypg
      external d01baw,d01bax,d01bay,d01baz
      
      if(itype.gt.4.or.itype.lt.1)then
       print*,'Getquad: itype is out of range'
       stop
      end if

      itypg=0.
  
      if(itype.eq.1)then
       call d01bbf(d01baz,a,b,itypg,n,weight,abscis,ifail)
      elseif(itype.eq.2)then
       call d01bbf(d01bay,a,b,itypg,n,weight,abscis,ifail)
      elseif(itype.eq.3)then
       call d01bbf(d01bax,a,b,itypg,n,weight,abscis,ifail)
      else
       call d01bbf(d01baw,a,b,itypg,n,weight,abscis,ifail)
      end if

      if(ifail.gt.0)then
       print*,'Getquad; ifail = ',ifail
       stop
      end if

      return

      end
