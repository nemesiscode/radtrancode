      real function goody_lor(s,w,v)
C     $Id: goody_lor.f,v 1.2 2011-06-17 15:40:26 irwin Exp $
C     *****************************************************************
C
C     Calculates the equivalent width of a band using Goody Model.
C 	Rodgers "Approx. Methods of Calculating Transmission' 1976
C
C     Assumes lines are Lorentz broadened
C   	Pat Irwin	18/2/94
C
C     *****************************************************************
      implicit none
      real s,w,v

      goody_lor = v/sqrt(1./(w*w) + 1./(s*s))

      return
      end
 
