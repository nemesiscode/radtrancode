      real function arctan(y,x)
C     **************************************************************
C     Function to find the unambiguous answer to tan-1(y/x).
C
C     Pat Irwin		1/1/05
C
C     **************************************************************
      real y,x,ang,pi

      parameter (pi=3.1415927)

      if(x.eq.0.0)then
       if(y.ge.0.0)then
        ang = 0.5*pi
       else
        ang = 1.5*pi
       endif 

       arctan = ang 

       return

      endif
  
      ang=atan(y/x)

      if(y.ge.0.)then
       if(x.ge.0.)then
        ang=ang
       else
        ang=ang+pi
       end if
      else
       if(x.ge.0.)then
        ang=2*pi+ang
       else
        ang=pi+ang
       end if
      end if

      arctan=ang

      return
      end
 
