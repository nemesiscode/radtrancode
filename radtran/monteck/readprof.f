      SUBROUTINE READPROF(IPFILE,NPRO,NGAS,RADIUS,MOLWT,XG,P,T,H,VMR,
     1  ID,ISO)
C     *************************************************************
C     Subroutine to read in a single atmospheric profile (.prf) file.
C
C     Input variables:
C	IPFILE		CHARACTER*100  Filename
C
C     Output variables:
C	NPRO		INTEGER	Number of points in profile.
C	NGAS		INTEGER	Number of gases.
C	RADIUS		REAL	Planetary radius (km)
C	MOLWT		REAL	Molecular weight (gm)
C	XG		REAL	Surface gravity (ms-2)
C	P(MAXPRO)		REAL	Pressure profile.
C	T(MAXPRO)		REAL	Temperature profile.
C	H(MAXPRO)		REAL	Height profile.
C	ID(MAXGAS)	INTEGER	Gas ID's
C	ISO(MAXGAS)	INTEGER	Gas ISO's
C
C     Pat Irwin		1/4/99
C  
C     *************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NPRO,NGAS,IFORM,N,I,J,K
      INTEGER IPLANET

      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL RADIUS,MOLWT,XGM,XCOEFF(3),LATITUDE
      REAL ELLIP,OMEGA,HEIGHT,XG
      INTEGER ID(MAXGAS),ISO(MAXGAS)
      CHARACTER*100 IPFILE,BUFFER
      CHARACTER*8 PNAME

C     ************************* CODE ******************************
      CALL REMSP(IPFILE)
      CALL FILE(IPFILE,IPFILE,'prf')
      CALL LOCASE(IPFILE)

      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
C     First skip header
1     FORMAT(A)
54    READ(1,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)IFORM
      READ(1,*)IPLANET,LATITUDE,NPRO,NGAS,MOLWT
      DO 20 I=1,NGAS
       READ(1,*)ID(I),ISO(I)
20    CONTINUE
C     reading the first block of profiles
      READ(1,*)
      N=MIN(NGAS,3)
C     N is the maximum VMR which can be read in from the next block
      DO 30 I=1,NPRO
       IF(IFORM.EQ.0)THEN
         READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,N)
       ELSE IF(IFORM.EQ.1)THEN
         READ(1,*)H(I),T(I),(VMR(I,J),J=1,N)
       ELSE   
         CALL WTEXT('invalid format')
         STOP   
       END IF
30    CONTINUE
C     reading in additional blocks if any
C     N VMR profiles have been read in so far
33    IF(NGAS.GT.N)THEN
        READ(1,*)
C       profiles up to VMR(?,K) to be read from this block
        K=MIN(NGAS,(N+6))
        DO 32 I=1,NPRO
         READ(1,*)(VMR(I,J),J=N+1,K)
32      CONTINUE
        N=K
        GOTO 33
      END IF

      CLOSE(UNIT=1)



      HEIGHT=0.0
      CALL NEWGRAV(IPLANET,LATITUDE,HEIGHT,RADIUS,XG,PNAME)

      PRINT*,'Planet is : ',PNAME
      PRINT*,'Radius, gravity is : ',RADIUS,XG

      RETURN
   
      END
