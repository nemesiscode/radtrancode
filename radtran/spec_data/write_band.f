      SUBROUTINE WRITE_BAND(POUT,WAVEN)
C     $Id: write_band.f,v 1.3 2011-06-17 15:53:03 irwin Exp $
C****************************************************************************
C_TITL:	WRITE_BAND
C
C_DESC:	Routine for writing a standard band parameter output file
C
C_ARGS:
C
C_CALL: prompt		Prompts the user for input.
C	file		Forces file extension.
C
C_HIST:	4/10/94	PGJI
C****************************************************************************

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/dbcom.f'

      INTEGER LUN
      PARAMETER (LUN=2)
      REAL POUT(MAXOUT,MAXGAS,7),WAVEN(MAXOUT,2)
      CHARACTER*100 BUFFER
      CHARACTER*78 LLINE,HEADER

      LLINE='----------------------------------------------------
     1--------------------------'

      HEADER='   V0       dV       Kv(T0)    delta/AD0      y0
     1           El       SFB      '

C********************************* CODE ********************************

      CALL PROMPT('Enter name of output file : ')
      READ(*,1)DBFILE
1     FORMAT(A)
      CALL FILE(DBFILE,DBFILE,'ban')

      WRITE(*,*)'%writing output'
      OPEN(UNIT=2,FILE=DBFILE,STATUS='UNKNOWN')
      WRITE(LUN,314)LLINE
      WRITE(LUN,302)DBFILE

      CALL PROMPT('Enter Originator : ')
      READ(*,1)BUFFER
      WRITE(LUN,307)BUFFER

      CALL PROMPT('Enter Date : ')
      READ(*,1)BUFFER
      WRITE(LUN,308)BUFFER

      BUFFER = ' '
      WRITE(LUN,303)BUFFER

      DO 305 I=1,NGAS
        WRITE(*,500)IDGAS(I),ISOGAS(I)
        READ(*,1)BUFFER
        WRITE(LUN,311)IDGAS(I),ISOGAS(I),BUFFER
305   CONTINUE
      PRINT*,'Enter any additional header information '
      READ(*,1)BUFFER
      WRITE(LUN,304)BUFFER
      BUFFER='*******************************************************'
      WRITE(LUN,304)BUFFER
      WRITE(LUN,401)VMIN
      WRITE(LUN,402)DELV
      WRITE(LUN,400)FWHM
      WRITE(LUN,403)NPOINT
      WRITE(LUN,404)NGAS
      DO 200 I=1,NGAS
        WRITE(LUN,405)IDGAS(I),ISOGAS(I)
200   CONTINUE

      WRITE(LUN,314)LLINE

      DO 205 I=1,NGAS
        WRITE(LUN,312)IDGAS(I),ISOGAS(I)
        WRITE(LUN,314)HEADER
        DO 210 J=1,NPOINT
          WRITE(LUN,932)WAVEN(J,1),WAVEN(J,2),POUT(J,I,1),POUT(J,I,2),
     1    POUT(J,I,3),POUT(J,I,4),POUT(J,I,5)
210     CONTINUE
205   CONTINUE

      CLOSE(LUN)

302   FORMAT(1X,'Original file name : ',1A56)
307   FORMAT(1X,'Originator : ',A60)
308   FORMAT(1X,'Date : ',A60)
303   FORMAT(1X,'Data sources : ',A60)
311   FORMAT(1X,'Gas : ',I3,I3,' ',A60)
304   FORMAT(1X,A60)
312   FORMAT(1X,'Band data for gas : ',I3,I3)
314   FORMAT(1X,A78)
400   FORMAT(1X,'FWHM = ',F8.2)
401   FORMAT(1X,'VMIN = ',F8.2)
402   FORMAT(1X,'DELV = ',F8.2)
403   FORMAT(1X,'NPOINT = ',I5)
404   FORMAT(1X,'NGAS = ',I3)
405   FORMAT(1X,'Gas : ',I3,I3)
500   FORMAT(1X,'Enter source of data for gas ',I3,I3)
932   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5)

      RETURN

      END
