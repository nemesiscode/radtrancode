      SUBROUTINE RDMOD(TEXT)
C     $Id:
C***********************************************************************
C_TITL:	RDMOD
C
C_DESC:	Reads in a model atmosphere for path.f
C
C_ARGS:	Input variable:
C	text	CHARACTER*(*)	Text string containing the name of the
C				profile runname.prf file. 
C				CHARACTER*(*) declares an incoming
C				character variable whose length is
C				unknown.

C_FILE:	unit=1	Atmospheric model (.prf)
C
C_CALL:	remsp		Remove spaces from input string.
C	file		Forces file extension.
C	locase		Make lower case the input string.
C	addgas		
C	newgrav		Sets the gravitational acceleration based on
C                       selected planet and latitude.
C
C_HIST:	26feb93	SBC	Original version.
C***********************************************************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C ../includes/pathcom.f holds the variables used by the software when
C calculating atmospheric paths (e.g. NPATH and IMOD).
      INCLUDE '../includes/laycom.f'
      INCLUDE '../includes/emcee.f'
C ../includes/laycom.f holds variables used only by the path software
C parameters are passed between routines mostly using common blocks
C because of the extensive use of large arrays. NOTE: laycom uses
C parameters defined in pathcom.

      INTEGER i,j,k,n,ilun
C N: Maximum VMR which can be read in from the next block.
C ILUN: File unit used for openning files.

      REAL g,dh,XFMIN,FRAC
C G: Gravitivational acceleration [m/s^2] determined as each height above
C the reference surface.
C DH: differential height between layers.
C XFMIN: How close the VMRs must add to 1.0 to pass inspection for AMFORM=1
      PARAMETER(XFMIN=0.001)
      CHARACTER*(*) text
      CHARACTER*8 pname
C PNAME: Planetary name.
      CHARACTER*100 buffer,ipfile
C IPFILE: Input filename.

C********************************* CODE ********************************

C Reading in in vertical profiles produced by profile.f

      ilun = 1
      READ(text,1)ipfile
1     FORMAT(A)
      CALL remsp(ipfile)

      CALL locase(ipfile)

      CALL file(ipfile,ipfile,'prf')
      print*,'len = ',len(ipfile)
      WRITE(*,*)' RDMOD.f :: reading model: ',ipfile
      OPEN(UNIT=ilun,FILE=ipfile,STATUS='OLD')
C First skip the header (if any)
54    READ(ILUN,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)AMFORM
      print*,'RDMOD: AMFORM = ',AMFORM
      IF(AMFORM.EQ.0)THEN
        READ(ILUN,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
      ELSE
        READ(ILUN,*)IPLANET,LATITUDE,NPRO,NVMR
      ENDIF


      IF(NVMR.GT.MAXGAS)THEN
       Print*,'Profile has too many gases. Either reduce or recompile'
       print*,'NVMR = ',NVMR
       print*,'MAXGAS = ',MAXGAS
      ENDIF

      DO 20 i=1,nvmr
        READ(ilun,*)id(i),iso(i)
C        print*,id(i),iso(i)
20    CONTINUE

      READ(ilun,*)
      DO 30 I=1,NPRO
        READ(ILUN,*)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
30    CONTINUE
      CLOSE(UNIT=82)

      IF(VMRflag.eq.1)THEN
       DO 120 I=1,NPRO
        DO 130 J=1, NVMR
         VMR(I,J) = MCMCvmr(J)
130     CONTINUE
        T(I) = MCMCtemp(i)
        H(I) = MCMCheight(i)
120    CONTINUE

      ENDIF

      MODEL = .TRUE.
      LAYERS = .FALSE.

      J=1
      DO 47 I=2,NPRO
        DH = H(I)-H(I-1)
        IF(DH.NE.0.0)THEN
          J=J+1
          H(J)=H(I)   
          P(J)=P(I)
          T(J)=T(I)
          DO K=1,NVMR   
            VMR(J,K)=VMR(I,K)
          ENDDO   
        ELSE
          WRITE(*,*)' RDMOD.f :: Identical heights found.'
          WRITE(*,*)' RDMOD.f :: Ignoring level: ',i
        ENDIF
47    CONTINUE

      IF(AMFORM.EQ.1)THEN
       DO 48 I=1,NPRO
        FRAC=0.
        DO J=1,NVMR
         FRAC=FRAC+VMR(I,J)
        END DO
        FRAC = ABS(FRAC-1.0)
        IF(FRAC.GT.XFMIN)THEN
         PRINT*,'Error in RDMOD.F. VMRs do not add up to 1.0'
        ENDIF
48     CONTINUE
      ENDIF

C Adding the model gases to the main arrays
      DO 21 I=1,NVMR
        CALL ADDGAS(ID(I),ISO(I),ATMGAS(I))
21    CONTINUE

      CALL NEWGRAV(IPLANET,LATITUDE,H(1),RADIUS,G,PNAME)

      RETURN

      END
