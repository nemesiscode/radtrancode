      double precision function get_albedo(nalb,valb,alb,x)
C     **********************************************************
C     Routine to interpolate a surface albedo spectrum to the required
C     wavenumber/wavelength.
C
C     Pat Irwin	17/11/05
C
C     **********************************************************
      integer nalb
      real valb(nalb),alb(nalb),x,y

      if(x.le.valb(1))then
       y = alb(1)
      else if(x.ge.valb(nalb))then
       y = alb(nalb)
      else
       call verint(valb,alb,nalb,y,x)
      endif

      get_albedo = dble(y)

      return

      end
