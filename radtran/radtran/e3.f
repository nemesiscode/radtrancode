      real function e3(x)
C     ****************************************************************
C     Utility routine to calculate the E3 exponential integral
C
C     Pat Irwin   30/9/02
C
C     ****************************************************************
      implicit none
      real x,e2

      if(x.eq.0)then
       e3 = 0.5
       return
      endif

      if(x.lt.8.0) then
        e3 = 0.5*(exp(-x) - x*e2(x)) 
      else
        e3 = exp(-x)*(1.0 -3.0/x + 3*4/x**2)/x
      endif

      return

      end

      real function e2(x)
C     ****************************************************************
C     Utility routine to calculate the E2 exponential integral
C
C     Pat Irwin   30/9/02
C
C     ****************************************************************
      implicit none
      real x,e1

      e2 = exp(-x)-x*e1(x)

      return

      end

      real function e1(x)
C     ****************************************************************
C     Utility routine to calculate the E1 exponential integral
C
C     Pat Irwin   30/9/02
C
C     ****************************************************************
      implicit none
      real cerr,sum,nfac,dsum,err,x,n,xfac
      cerr = 0.0001             ! % convergence limit

      sum = -0.5772156 -log(abs(x))

      n=0.0
      xfac = -1.0
10    n = n + 1.0
      xfac = -xfac*(x/n)
      dsum = xfac/n
      sum=sum+dsum

      err = 100*abs(dsum/sum)

      if(err.gt.cerr) goto 10

      e1 = sum

      return

      end   
