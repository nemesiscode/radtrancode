************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C                          SUBROUTINE INIT_SOLAR_WAVE
C
C
C	Calculates solar/stellar flux using a look up table. The table is read
C	only once per loop through the bins, with values stored in the
C	common block.
C
C       Input variables
C	ispace	integer	Wavespace of calculation
C	aname	character*100	Filename to read
C
C       The required units of the spectrum read in are:
C        Wavenumber:	W (cm-1)-1, i.e. the total spectral power output
C        Wavelength:	W um-1 i.e. the total spectral power output
C	Both these quantities should be calculated as I*4*PI*R^2 where
C	 I is the surface spectral irradiance of the sun/star (i.e W m-2 
C           (cm-1)-1 or W m-2 um-1)
C	 R is the radius of the Sun/Star (units of m)
C
C	Pat Irwin	17/6/11	Original Version
C	Pat Irwin	1/3/12	Updated for Radtrans2.0
C
C-----------------------------------------------------------------------

	SUBROUTINE init_solar_wave (ispace,aname)

	IMPLICIT NONE
        include '../includes/arrdef.f'

	INTEGER		iunit, npt, I,ispace,ispace1,iread
	REAL		wave(maxbin), rad(maxbin), y
        CHARACTER*100	aname,solfile
        CHARACTER*80    dummy
        PARAMETER (iunit=26)

        common/solardat/iread,wave, rad,  npt

C        call datarchive(aname)
        solfile = aname

1       format(a)

	open (iunit,file=solfile,status='old')

C       Skip header
54      read(iunit,1)dummy
        if(dummy(1:1).eq.'#')goto 54

        read(dummy,*)ispace1

        if(ispace1.ne.ispace)then
         print*,'Error in INIT_SOLAR_WAVE'
         print*,'wavespace of solar/stellar file is not the same'
         print*,'as the calculation wavespace'
         print*,ispace,ispace1
         stop
        endif

	npt = 0
        do i = 1,maxbin
  		  read (iunit,*,end=20) wave(i), rad(i)
		  npt = npt + 1
        enddo
20	continue
  	close (iunit)

        iread=999

        return

        end
