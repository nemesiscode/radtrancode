      SUBROUTINE READPRFHEIGHT(IPFILE,NPRO,HEIGHT)
C     $Id:
C     ***********************************************************************
C     Subroutine to read a .prf file and return the number of vertical
C     levels and the height array.
C
C     Input filename
C	IPFILE	CHARACTER*100 	Input filename
C
C     Output filename
C    	NPRO	INTEGER		Number of levels
C	HEIGHT(NPRO)	REAL	Heights
C
C     Pat Irwin 17/10/03	Original
C     Pat Irwin 4/6/14		Modified from readrefhead.f
C
C     ***********************************************************************
      IMPLICIT NONE
      include '../includes/arrdef.f'

      CHARACTER*100 IPFILE,BUFFER
      CHARACTER*8 PNAME
      INTEGER NPRO,NVMR,AMFORM,IPLANET,ID(MAXGAS),ISO(MAXGAS),I
      REAL HEIGHT(MAXPRO),P(MAXPRO),T(MAXPRO)
      REAL LATITUDE,MOLWT
      LOGICAL GASGIANT

      CALL FILE(IPFILE,IPFILE,'prf')
      
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
C     First skip header
54     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 54
       READ(BUFFER,*)AMFORM
1      FORMAT(A)
       IF(AMFORM.EQ.0)THEN
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       ELSE
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
       ENDIF
       DO 20 I=1,NVMR
         READ(1,*)ID(I),ISO(I)
20     CONTINUE
C      Skip header
       READ(1,*)
       DO 30 I=1,NPRO
         READ(1,*)HEIGHT(I),P(I),T(I)
30     CONTINUE

      CLOSE(UNIT=1)

      RETURN

      END

