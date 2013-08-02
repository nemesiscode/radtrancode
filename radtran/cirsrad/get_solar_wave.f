************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C                          SUBROUTINE GET_SOLAR_WAVE
C
C
C	Calculates solar flux using a look up table. The table is read
C	only once per loop through the bins, with values stored in the
C	common block.
C
C-----------------------------------------------------------------------

	SUBROUTINE get_solar_wave (x, dist, solar)

	IMPLICIT NONE
        include '../includes/arrdef.f'
        include '../includes/constdef.f'

	INTEGER		npt,iread,iform
	REAL		x, dist, solar, wave(maxbin), rad(maxbin), y
        REAL		area, solrad

        common/solardat/iread, iform, solrad, wave, rad, npt

        if(iread.ne.999)then
         print*,'GET_SOLAR_WAVE: error. No solar spectrum initialised'
         stop
        endif

	call interp(wave, rad, npt, y, x)

C       If dist is > 0 then we calculate the flux at the distance of the 
C       planet from the star. Hence we need to divide by 4*pi*dist^2
C       to provide solar flux in units of W cm-2 um-1 or W cm-2 (cm-1)-1.
C      
C       If dist is < 0, then we want the total spectral irradiance
C       in units of W um-1 or W (cm-1)-1
       
        if(dist.gt.0.)then
          area = 4*PI*(dist*AU)**2
          solar = y/area
        else
          solar = y
        endif

	return

	end

************************************************************************
************************************************************************


