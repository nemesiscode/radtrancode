************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C_TITLE:		SUBROUTINE CIRSrtf
C
C_ARGS: Input Variables
C
C	CONVOUT
C       DIST:REAL       Distance from Sun (as prompted by CIRSDRV) in
C                       units of AU.
C       I:INT           Incrementor/Loop Counter for
C       INORMAL:INT     flag for ortho:para ratio (0=equilibrium =>1:1)
C                       (1=normal =>3:1)
C       ITYPE:INT       Value designating the chosen scattering routine
C                       (currently only scloud8 through scloud11).
C       J:INT           Incrementor/Loop counter for
C       K:INT           Incrementor/Loop counter for
C       NCONV:INT
C       NWAVE:INT       The number of calculation wavenumbers.
C       NPATH:INT       Number of individual paths through the layers.
C       OPFILE:CHARA*100	Operation filename
C	VCONV:REAL
C	VWAVE:REAL      Bin centres in wavenumber space.
C	Y
C	YOUT
C
C
C       Output variables
C       
C       OUTPUT          Output values at each wavenumber for each output
C                       type for each path
C
C_DESC:
C
C_HIST: 
C-----------------------------------------------------------------------

        SUBROUTINE cirsrtf(opfile, Dist, INormal, FWHM, vwave, nwave,
     1          npath, output, vconv, nconv, itype, convout)

	IMPLICIT NONE

C       Defines the maximum values for a series of variables (layers,
C         bins, paths, etc.)
	INCLUDE '../includes/arrdef.f'

        CHARACTER*100	opfile,FWHMFILE
        INTEGER         nwave, nconv, npath, itype, I, J, K
	INTEGER		INormal
        INTEGER         NFWHM,MFWHM
        PARAMETER       (MFWHM=1000)
	REAL		Dist, FWHM
        REAL          vwave(nwave), vconv(nconv), convout(maxout3),
     1                  output(maxout3), y(maxout), yout(maxout)
        LOGICAL         FWHMEXIST
        REAL            VFWHM(MFWHM),XFWHM(MFWHM)

       
        PRINT*,'CIRSRTF Called'
C       See if file is present forcing FWHM to vary with wavelength/wavenumber
        CALL FILE(OPFILE,FWHMFILE,'fwh')
        INQUIRE(FILE=FWHMFILE,EXIST=FWHMEXIST)
C       If such a file exists then read in the data
        IF(FWHMEXIST)THEN
         print*,'Reading FWHM information from : ',FWHMFILE
         OPEN(13,FILE=FWHMFILE,status='old')
          READ(13,*)NFWHM
          DO I=1,NFWHM
           READ(13,*)VFWHM(I),XFWHM(I)
          ENDDO
         CLOSE(13)
        ENDIF

C-----------------------------------------------------------------------
C
C       Call subroutine subCIRSrtf:
C
C-----------------------------------------------------------------------



        CALL subcirsrtf(opfile, Dist, INormal, vwave, nwave, npath,
     1		itype, output)

C-----------------------------------------------------------------------
C       We need to check that MAXOUTPUT is big enough to allow CONVOUT   
C       to contain all of the convolved output, ie we must check that   
C                 NCONV * NPATH  <=  MAXOUT3
C       It would be better to do this before the calculation, but NPATH
C       is first determined in subnrtf, and subnrtf does not know the
C       value of NCONV.
C-----------------------------------------------------------------------

        IF ((NCONV * NPATH).GT.MAXOUT3) THEN
           PRINT*, 'FATAL ERROR: NCONV*NPATH > MAXOUT3, so that the '
           PRINT*,'convolved output will not fit into the output array.'
           PRINT*, 'NCONV= ', nconv
           PRINT*, 'NPATH= ', npath
           PRINT*, 'MAXOUT3= ', maxout3
           STOP
        ENDIF

        DO I= 1, npath
           DO J= 1, nwave
              K= I+(J-1)*npath
              y(J)= output(K)
           ENDDO


           CALL cirsconv(opfile,fwhm, nwave,vwave,y,nconv,vconv,yout,
     1 FWHMEXIST,NFWHM,VFWHM,XFWHM)

           DO J= 1, nconv
              K= I+(J-1)*npath
              convout(K)= yout(J)
           ENDDO

        ENDDO


	RETURN

	END

************************************************************************
************************************************************************
