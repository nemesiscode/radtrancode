************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C_TITLE:		SUBROUTINE CIRSrtf_waveMC
C
C_ARGS: Input Variables
C
C	CONVOUT
C       DIST:REAL       Distance from Sun (as prompted by CIRSDRV) in
C                       units of AU.
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1) 
C       I:INT           Incrementor/Loop Counter for
C       INORMAL:INT     flag for ortho:para ratio (0=equilibrium =>1:1)
C                       (1=normal =>3:1)
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

        SUBROUTINE cirsrtf_waveMC(opfile1, Dist, INormal, Iray, FWHM, 
     1          ispace, vwave, nwave, sol_ang,emiss_ang,aphi,
     2          npath, output, vconv, nconv, nem, vem, 
     3          emissivity,tsurf, convout)

	IMPLICIT NONE

C       Defines the maximum values for a series of variables (layers,
C         bins, paths, etc.)
	INCLUDE '../includes/arrdef.f'

        CHARACTER*100	opfile,opfile1,FWHMFILE
        INTEGER         nwave, nconv, npath, I, J, K
	INTEGER		INormal,Iray, ispace,nem,NFWHM,MFWHM
        PARAMETER	(MFWHM=1000)
	REAL		Dist, FWHM,sol_ang,emiss_ang,aphi
        REAL          vwave(nwave), vconv(nconv), convout(maxout3),
     1                  output(maxout3), y(maxout), yout(maxout),tsurf,
     2			vem(MAXSEC),emissivity(MAXSEC)
        REAL		VFWHM(MFWHM),XFWHM(MFWHM)
        LOGICAL		FWHMEXIST
C-----------------------------------------------------------------------
C
C       Call subroutine subCIRSrtf_wave:
C
C-----------------------------------------------------------------------

        opfile=opfile1

C       See if file is present forcing FWHM to vary with wavelength/wavenumber
        CALL FILE(OPFILE,FWHMFILE,'fwhm')
        INQUIRE(FILE=FWHMFILE,EXIST=FWHMEXIST)
C       If such a file exists then read in the data
        IF(FWHMEXIST)THEN
         OPEN(13,FILE=FWHMFILE,status='old')
          READ(13,*)NFWHM
          DO I=1,NFWHM
           READ(13,*)VFWHM(I),XFWHM(I)
          ENDDO
         CLOSE(13)
        ENDIF

        print*,'cirsrtf, opfile = ',opfile
        CALL subcirsrtf_waveMC(opfile, Dist, INormal, Iray, ispace, 
     1          vwave,nwave,npath,sol_ang, emiss_ang,aphi,nem,vem,
     2          emissivity,tsurf,output)

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


           CALL cirsconv(opfile,fwhm, nwave, vwave, y, nconv, 
     1      vconv,yout,FWHMEXIST,NFWHM,VFWHM,XFWHM)

           DO J= 1, nconv
              K= I+(J-1)*npath
              convout(K)= yout(J)
           ENDDO

        ENDDO


	RETURN

	END

************************************************************************
************************************************************************
