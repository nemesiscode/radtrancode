      SUBROUTINE READREFIPLAN(IPFILE,IPLANET,XLAT,RADIUS)
C     $Id:
C     ***********************************************************************
C     Subroutine to read a .ref file and return the number of vertical
C     levels and number of gases contained.
C
C     Input filename
C	IPFILE	CHARACTER*100 	Input filename
C	XLAT	REAL		Required latitude
C
C     Output filename
C    	IPLANET	INTEGER		Planet identifier
C	RADIUS	REAL		Radius (km) at 0km altitude level
C
C     Pat Irwin 21/07/10	Modified from readrefhead.pro
C
C     ***********************************************************************
      IMPLICIT NONE
      CHARACTER*100 IPFILE,BUFFER
      INTEGER NPRO,NVMR,AMFORM,IPLANET,NLATREF
      REAL LATITUDE,MOLWT,RADIUS,G,H1,XLAT
      LOGICAL GASGIANT
      CHARACTER*8 PNAME

      CALL FILE(IPFILE,IPFILE,'ref')

      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
C     First skip header
54     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 54
       READ(BUFFER,*)AMFORM
1      FORMAT(A)
       READ(1,*)NLATREF
       IF(AMFORM.EQ.0)THEN
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       ELSE
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
       ENDIF
      CLOSE(UNIT=1)


      H1=0.0
      CALL NEWGRAV(IPLANET,XLAT,H1,RADIUS,G,PNAME)
      RADIUS=RADIUS

      RETURN

      END

