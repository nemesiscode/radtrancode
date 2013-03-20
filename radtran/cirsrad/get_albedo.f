      double precision function get_albedo(nem,vem,emissivity,x)
C     **********************************************************
C     Routine to interpolate a surface albedo spectrum to the required
C     wavenumber/wavelength.
C
C     Pat Irwin	17/11/05 Original
C     Pat Irwin 20/3/13  Revised to get albedo from emissivity spectrum
C
C     **********************************************************
      integer nem
      real vem(nem),emissivity(nem),x,y

      if(x.le.vem(1))then
       y = emissivity(1)
      else if(x.ge.vem(nem))then
       y = emissivity(nem)
      else
       call verint(vem,emissivity,nem,y,x)
      endif

      get_albedo = dble(1.0-y)

      return

      end
