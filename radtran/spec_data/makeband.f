      PROGRAM MAKEBAND
C     $Id: makeband.f,v 1.3 2011-06-17 15:53:01 irwin Exp $
C***********************************************************************
C_TITL:	MAKEBAND
C
C_DESC:	Reads in the intervals and number of gases to start the subroutine
C	LBLBAND. LBLBAND reads in the linedata and calculates the 
C	statistical random band parameters for each interval.
C
C_ARGS:
C
C_FILE:	UNIT LUN=2
C
C_CALL:	RDKEY
C	RDGAS
C	RDISO
C	LBLBAND
C	WRITE_BAND
C
C_HIST:	24/11/93	PGJI	Adapted from the GENLBL program LBL.
C	1/6/94		PGJI	Adapted from GENBAND
C	4/10/94		PGJI	Adapted from LBL_PARAMETER
C***********************************************************************

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/dbcom.f'

      INTEGER LUN
      PARAMETER (LUN=2)
C MAXOUT the maximum number of output bins
      REAL POUT(MAXOUT,MAXGAS,7),WAVEN(MAXOUT,2)

C********************************* CODE ********************************

      CALL PROMPT('Calculate from linedata (1) or set to zero (2)? ')
      READ*,ISPEC

      IF(ISPEC.EQ.1)THEN

        CALL PROMPT('Enter name of lbl key file : ')
        READ(5,1)KEYFIL
        CALL FILE(KEYFIL,KEYFIL,'key')

        CALL RDKEY(LUN)
        CALL RDGAS
        CALL RDISO
      ENDIF

      PRINT*,'Enter wavenumber range and spacing.'
      CALL PROMPT('Enter Vmin, Vmax, DelV, FWHM : ')
      READ*,VMIN,VMAX,DELV,FWHM
      NPOINT=1+INT((VMAX-VMIN)/DELV)

      print*,vmin,vmax,delv,fwhm,npoint

      IF(NPOINT.GT.MAXOUT)THEN
        WRITE(*,541)
541     FORMAT(' NPOINT>MAXOUT Too many bins. Recompile')
        STOP
      ENDIF

      CALL PROMPT('Enter number of gases : ')
      READ*,NGAS

      IF(NGAS.GT.MAXGAS)THEN
        WRITE(*,542)
542     FORMAT(' NGAS>MAXGAS Too many gases. Recompile')
      END IF

      DO J=1,NGAS
        WRITE(*,*)'Gas : ',J
        CALL PROMPT('Enter Gas ID and ISO : ')
        READ*,IDGAS(J),ISOGAS(J)
      ENDDO

      IF(ISPEC.EQ.1)THEN
        CALL LBLBAND(POUT,WAVEN)
      ELSE
        DO 733 I=1,NPOINT
          WAVEN(I,1)=VMIN+(I-1)*DELV
          WAVEN(I,2)=FWHM
          DO 730 IGAS=1,NGAS
            POUT(I,IGAS,1)=0.
            DO 732 J=2,5
              POUT(I,IGAS,J)=1.0
732         CONTINUE
730       CONTINUE
733     CONTINUE
      ENDIF
      CALL WRITE_BAND(POUT,WAVEN)

      CALL FILE(DBFILE,KEYFIL,'key')

      OPEN(12,FILE=KEYFIL,STATUS='UNKNOWN')
      WRITE(12,1)DBFILE
      WRITE(12,1)GASFIL
1     FORMAT(A)
      CLOSE(12)

      WRITE(*,*)' MAKEBAND.f :: calculation complete.'

      END
