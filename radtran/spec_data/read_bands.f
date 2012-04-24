      SUBROUTINE READ_BANDS
C     $Id: read_bands.f,v 1.6 2011-06-23 09:08:15 irwin Exp $
C     **********************************************************************
C     Routine to read in a band parameter file.
C
C     Pat Irwin		4/10/94
C     **********************************************************************
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'
C     **********************************************************************
      INTEGER LUN,I,OFFSET
      PARAMETER (LUN=2)
      REAL WAVEN(MAXBIN+300,2)
      REAL TPOUT(MAXBIN+300,MAXGAS,7),TWAVEN(MAXBIN+300,2)

      CHARACTER*20 HEAD
      CHARACTER*10 BUFFER
      CHARACTER*100 LINE,LINE1
      REAL VMIN1,DELV1,VMA,VMA1,FWHM1,READL
      INTEGER NPOINT1,NGAS1,GIDENT(MAXGAS),IBAND(MAXGAS)
      INTEGER IDGAS1(MAXGAS),ISOGAS1(MAXGAS)
      LOGICAL FLAG
C     **********************************************************************

      CALL FILE(DBFILE,DBFILE,'ban')
      
      OPEN(UNIT=2,FILE=DBFILE,STATUS='OLD')

C     FIRST SKIP HEADER
1     FORMAT(A)
11    READ(LUN,500)BUFFER
      IF(BUFFER.NE.'**********')GOTO 11

      READ(LUN,1)LINE
      VMIN1 = READL(LINE)
      READ(LUN,1)LINE
      DELV1 = READL(LINE)
      READ(LUN,1)LINE
      FWHM1 = READL(LINE)
      READ(LUN,1)LINE
      NPOINT1 = INT(READL(LINE))
      READ(LUN,1)LINE
      NGAS1 = INT(READL(LINE))

C      READ(LUN,401)VMIN1
C      READ(LUN,402)DELV1
C      READ(LUN,400)FWHM1
C      READ(LUN,403)NPOINT1
C      READ(LUN,404)NGAS1

      print*,vmin,delv,fwhm,npoint,ngas
      print*,vmin1,delv1,fwhm1,npoint1,ngas1

      VMA=VMIN+(NPOINT-1)*DELV
      VMA1=VMIN1+(NPOINT1-1)*DELV1

      IF(DELV.NE.DELV1)THEN
       PRINT*,'Read_bands: Wave step DELV incompatible'
       PRINT*,'DELV, DELV1 = ',DELV,DELV1
       STOP
      END IF

      IF(FWHM.NE.FWHM1)THEN
       PRINT*,'Read_bands: Resolution FWHM incompatible'
       PRINT*,'FWHM, FWHM1 = ',FWHM,FWHM1
       STOP
      END IF

      IF(VMIN.LT.VMIN1.OR.VMA.GT.VMA1)THEN
       PRINT*,'Read_bands: Band data does not cover entire range'
       PRINT*,'VMIN, VMAX = ',VMIN,VMA
       PRINT*,'VMIN_TABLE, VMAX_TABLE = ',VMIN1,VMA1
      END IF

      IF(NGAS.GT.NGAS1)THEN
       PRINT*,'Read_bands: Warning, not enough gases in band file'
       PRINT*,'NGAS, NGAS1 = ',NGAS,NGAS1
      END IF

      IF(NGAS.GT.MAXBGAS)THEN
       PRINT*,'Read_bands: NGAS > MAXBGAS'
       PRINT*,NGAS,MAXBGAS
       PRINT*,'Reduce NGAS, or increase MAXBGAS and recompile'
       STOP
      ENDIF 

      DO 200 I=1,NGAS1
         READ(LUN,1)LINE
         J=LEN(LINE)
         LINE1 = LINE(7:J)
         READ(LINE1,*)IDGAS1(I),ISOGAS1(I)

C        READ(LUN,405)IDGAS1(I),ISOGAS1(I)
C        print*,idgas1(i),isogas1(i)
200   CONTINUE

      DO 201 I=1,NGAS
       FLAG=.FALSE.
       DO 202 J=1,NGAS1
        IF(IDGAS(I).EQ.IDGAS1(J).AND.ISOGAS(I).EQ.ISOGAS1(J))THEN
         FLAG=.TRUE.
         GIDENT(I)=J
        END IF
202    CONTINUE
       IF(.NOT.FLAG)THEN
        PRINT*,'read_band : Warning, Gas ',I,' not in file'
        GIDENT(I)=0
       END IF
       print*,'id iso gident: ',idgas(i),isogas(i),gident(i)       
201   CONTINUE   

C     Read in end of header

      READ(LUN,406)HEAD

      DO 205 I=1,NGAS1
       READ(LUN,406)HEAD
       READ(LUN,406)HEAD
C       print*,HEAD
C       print*,HEAD(1:3)
       IF(HEAD(1:3).EQ.'Ex.')THEN
        DO 211 J=1,NPOINT1
         READ(LUN,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         READ(LUN,933)TPOUT(J,I,6),TPOUT(J,I,7)
211     CONTINUE
        IBAND(I)=1
       ELSEIF(HEAD(1:3).EQ.'2-E')THEN
        DO 209 J=1,NPOINT1
         READ(LUN,934)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5),TPOUT(J,I,6),TPOUT(J,I,7)
209     CONTINUE
        IF(HEAD(1:4).EQ.'2-E1')THEN
         IBAND(I)=2
        ELSE
         IBAND(I)=3
        ENDIF
       ELSEIF(HEAD(1:4).EQ.'Kark')THEN
        DO 179 J=1,NPOINT1
         READ(12,*)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4)
179     CONTINUE
        IBAND(I)=4
       ELSE
        DO 210 J=1,NPOINT1
         READ(LUN,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         TPOUT(J,I,6)=0.0
         TPOUT(J,I,7)=0.0
         IBAND(I)=0
210     CONTINUE
       ENDIF
205   CONTINUE

      CLOSE(LUN)

      OFFSET = INT((VMIN - VMIN1)/DELV)
      print*,'Data read. offset = ',offset

      DO 206 J=1,NGAS
       IF(GIDENT(J).GT.0)THEN
        BANDTYP(J)=IBAND(GIDENT(J))
       ELSE
        BANDTYP(J)=0
       ENDIF
       print*,'Gas ID, Band type : ',J,BANDTYP(J)
      
       DO 207 I=1,NPOINT
        DO 208 K=1,7
         BANDPAR(I,J,K)=0.
         IF(GIDENT(J).GT.0)THEN
          IF((I+OFFSET).GT.0.0.AND.(I+OFFSET).LE.NPOINT1)THEN
           BANDPAR(I,J,K)=TPOUT(I+OFFSET,GIDENT(J),K)
          ENDIF
         END IF
208     CONTINUE
207    CONTINUE
206   CONTINUE

500   FORMAT(1X,A10)
401   FORMAT(1X,'VMIN = ',F8.2)
402   FORMAT(1X,'DELV = ',F8.2)
400   FORMAT(1X,'FWHM = ',F8.2)
403   FORMAT(1X,'NPOINT = ',I5)
404   FORMAT(1X,'NGAS = ',I3)      
405   FORMAT(1X,'Gas : ',I3,I3)
406   FORMAT(1X,A)
932   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5)
933   FORMAT(1x,2(e12.5))
934   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5)

      RETURN

      END



      REAL FUNCTION READL(LINE)
      CHARACTER*100 LINE,STRIP
      INTEGER J,I
      REAL X
      
      
      DO I=1,10
       IF(LINE(I:I).EQ.'=')J=I
      ENDDO
      STRIP = LINE(J+1:30)

      READ(STRIP,*)X

      READL=X

      RETURN

      END
