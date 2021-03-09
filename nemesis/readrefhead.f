      SUBROUTINE READREFHEAD(IPFILE,NPRO,NVMR,GASGIANT)
C     $Id:
C     ***********************************************************************
C     Subroutine to read a .ref file and return the number of vertical
C     levels and number of gases contained.
C
C     Input filename
C	IPFILE	CHARACTER*100 	Input filename
C
C     Output filename
C    	NPRO	INTEGER		Number of levels
C	NVMR	INTEGER		Number of gases
C	GASGIANT LOGICAL	Is the planet a Gas Giant?
C
C     Pat Irwin 17/10/03	Original
C
C     ***********************************************************************
      IMPLICIT NONE
      CHARACTER*100 IPFILE,BUFFER
      CHARACTER*8 PNAME
      INTEGER NPRO,NVMR,AMFORM,IPLANET,ISURF,NLATREF
      REAL LATITUDE,MOLWT
      LOGICAL GASGIANT
      REAL XGM,XCOEFF(3),XRADIUS,XELLIP,XOMEGA
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


      CALL FILE(IPFILE,IPFILE,'ref')
      
      if(idiag.gt.0)print*,'Reading reference file from : ',IPFILE
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

      GASGIANT=.TRUE.

      CALL READ_GRAV(IPLANET,XGM,XCOEFF,XRADIUS,XELLIP,XOMEGA,PNAME,
     1 ISURF)

      IF(ISURF.EQ.1)GASGIANT=.FALSE.
      
      RETURN

      END

