c-----------------------------------------------------------------------
      SUBROUTINE zgauleg(x,w,n,np)
c-----------------------------------------------------------------------
c
c	n = number of points
C	np = array size
c	x = g_ord
c	w = del_g
c
c-----------------------------------------------------------------------
c	4/4/06	nick teanby	original code (based on num rec gauleg)
c-----------------------------------------------------------------------
      INTEGER n
      REAL x1,x2,x(np),w(np)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      
c-----------------------------------------------------------------------
      x1=0.
      x2=1.
	do i=1,np
        x(i) = 0.
        w(i) = 0.
      enddo
c-----------------------------------------------------------------------
      
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.
