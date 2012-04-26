      PROGRAM NEMESISPROFILE
C     $Id%
C--------------------------------------------------------------
      IMPLICIT NONE
      INTEGER MAXPRO,MAXVMR
      PARAMETER(MAXPRO=300,MAXVMR=20)
C     MAXPRO is the maximum number of vertical points which can be stored.
C     MAXVMR is the maximum number of vertical mixing ratio profiles.
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXVMR),TINY
C     TINY is the smallest allowable VMR (NEMESIS uses log amounts)
      PARAMETER(TINY=1E-32)
      INTEGER NPRO,NVMR,ID(MAXVMR),ISO(MAXVMR)
C     H is the height in kilometres above some NOMINAL zero.
C     P is the pressure in atmospheres (not bar).
C     T is the temperature in Kelvin.
C     VMR holds the NVMR volume mixing ratio profiles for each of the gases.
C     There are NPRO points in each profile.
C     ID and ISO hold the local identifier and isotope identifier for
C     each gas. Note that this program does not check that you only include
C     each gas once or that the identifiers are valid.
C
C---------------------------------------------------------------------------- 
      CHARACTER*100 IPFILE,OPFILE,BUFFER
      INTEGER I,J,K,N,IFORM,IPLANET
      INTEGER IPROF,IGAS,IERR
      REAL TEMP,XERR,CLEN
      REAL LATITUDE,X1,E1
      REAL MOLWT
C----------------------------------------------------------------------------
C
C     First reading in an existing set of profiles which could have been
C     produced using a text editor or a previous use of this program or a
C     project specific program.
C     The nominal format (IFORM=0) is:
C     IFORM
C     IPLANET,LATITUDE,NPRO NVMR MOLWT
C     ID(1) ISO(1)
C     :
C     :
C     ID(NVMR) ISO(NVMR)
C     "H        P           T           gas1         gas2         gas3"
C      H(1)     P(1)        T(1)        VMR(1,1)     VMR(1,2)     VMR(1,3)
C      :        :           :           :            :            :
C      :        :           :           :            :            :
C      H(NPRO   P(NPRO)     T(NPRO)     VMR(NPRO,1)  VMR(NPRO,2)  VMR(NPRO,3)
C     "gas4        gas5.....     gasVMR   "
C      VMR(1,4)    VMR(1,5)......VMR(1,NVMR)
C      :           :             :
C      :           :             : 
C      VMR(NPRO,4) VMR(NPRO,5)   VMR(NPRO,NVMR)
C     i.e. profiles are stored in blocks of 6 columns each with a descriptive
C     header. Other value of IFORM allow minor variations to format, mainly
C     for data input.
C     IFORM=1 assumes pressure profile unknown
       CALL PROMPT('Enter filename : ')
       READ(*,10)IPFILE
10     FORMAT(A)
       CALL FILE(IPFILE,IPFILE,'prf')
       OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
C     First skip header
54     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 54
       READ(BUFFER,*)IFORM
1      FORMAT(A)
       READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       DO 20 I=1,NVMR
       READ(1,*)ID(I),ISO(I)
20     CONTINUE
C      reading the first block of profiles
       READ(1,*)
       N=MIN(NVMR,3)
C      N is the maximum VMR which can be read in from the next block
       DO 30 I=1,NPRO
       IF(IFORM.EQ.0)THEN
         READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,N)
        ELSE IF(IFORM.EQ.1)THEN
         READ(1,*)H(I),T(I),(VMR(I,J),J=1,N)
        ELSE
         CALL WTEXT('invalid format')
         STOP
         END IF
30     CONTINUE
C      reading in additional blocks if any
C      N VMR profiles have been read in so far
33     IF(NVMR.GT.N)THEN 
	READ(1,*)
C       profiles up to VMR(?,K) to be read from this block
	K=MIN(NVMR,(N+6))
	DO 32 I=1,NPRO
	READ(1,*)(VMR(I,J),J=N+1,K)
32      CONTINUE
	N=K
	GOTO 33
       END IF
       CLOSE(UNIT=1)

C      all processing below assumes that heights are in ascending order
C      so sorting just in case
       DO 12 J=1,NPRO
       DO 12 I=1,NPRO-1
       IF(ABS(H(I)-H(I+1)).LT.0.01)THEN
	WRITE(*,14)
14      FORMAT(' identical height values found')
C	STOP
       END IF
       IF(H(I).GT.H(I+1))THEN
	TEMP=H(I+1)
	H(I+1)=H(I)
	H(I)=TEMP
	TEMP=P(I+1)
	P(I+1)=P(I)
	P(I)=TEMP
	TEMP=T(I+1)
	T(I+1)=T(I)
	T(I)=TEMP
	DO 15 K=1,NVMR
 	 TEMP=VMR(I+1,K)
	 VMR(I+1,K)=VMR(I,K)
	 VMR(I,K)=TEMP
15      CONTINUE
       END IF
12     CONTINUE

      CALL PROMPT('Extract T(0) or vmr (1) profile : ')
      READ*,IPROF

      IF(IPROF.EQ.1)THEN
       CALL PROMPT('Enter IGAS : ')
       READ*,IGAS
      ENDIF

      CALL PROMPT('Constant error (0) or percentage error (1) : ')
      READ*,IERR
      
      CALL PROMPT('Enter error value : ')
      READ*,XERR

      CALL PROMPT('Enter correlation length (logp units) : ')
      READ*,CLEN


     
      CALL PROMPT('Enter output filename : ')
      READ(5,1)OPFILE

      OPEN(13,FILE=OPFILE,STATUS='UNKNOWN')
       WRITE(13,*)NPRO,CLEN 
       DO 22 I=1,NPRO 
        IF(IPROF.EQ.0)THEN
         X1 = T(I)
        ELSE
         X1 = VMR(I,IGAS)
        ENDIF
C       Make sure value is not too small
        IF(X1.LT.TINY)X1=TINY
        IF(IERR.EQ.0)THEN
         E1 = XERR
        ELSE
         E1 = XERR*X1/100.0
        ENDIF

        WRITE(13,*)P(I),X1,E1
     
22     CONTINUE

      CLOSE(13)

      END
