      real function malkmus_lor(s,w,v)
C     $Id: malkmus_lor.f,v 1.2 2011-06-17 15:40:27 irwin Exp $
C     *****************************************************************
C
C     Calculates the equivalent width of a band using Malkmus Model.
C       Rodgers "Approx. Methods of Calculating Transmission' 1976
C
C     Assumes lines are Lorentz broadened
C       Pat Irwin       18/2/94
C
C     *****************************************************************

      implicit none
      real s,w,a,v


      a=1. + 4.*w*w/(s*s)
      a=sqrt(a) - 1.
      malkmus_lor = v*s*s*a*0.5/w


      return
      end
 
