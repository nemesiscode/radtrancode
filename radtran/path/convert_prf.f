      PROGRAM CONVERT_PRF
C     *********************************************************************
C     Program to convert old-style 'block' .prf files to new format with all 
C     variables for a given level on the same row.
C
C     Pat Irwin		1/4/14
C
C     *********************************************************************
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'

C     MAXPRO is the maximum number of vertical points which can be stored.
C     MAXGAS is the maximum number of vertical mixing ratio profiles.
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      INTEGER NPRO,NVMR,ID(MAXGAS),ISO(MAXGAS)
C     H is the height in kilometres above some NOMINAL zero.
C     P is the pressure in atmospheres (not bar).
C     T is the temperature in Kelvin.
C     VMR holds the NVMR volume mixing ratio profiles for each of the gases.
C     There are NPRO points in each profile.
C     ID and ISO hold the local identifier and isotope identifier for
C     each gas. Note that this program does not check that you only include
C     each gas once or that the identifiers are valid.
      INTEGER IPLANET,AMFORM
      REAL MOLWT
      CHARACTER*100 IPFILE,OPFILE,BUFFER

C     First reading in an existing set of profiles which could have been
C     produced using a text editor or a previous use of this program or a
C     project specific program.
C     The nominal format (AMFORM=0) is:
C     AMFORM
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
C     header. Other value of AMFORM allow minor variations to format, mainly
C     for data input.
C
C     When AMFORM=1, it is assumed that the vmrs add up to 1.0 and this the
C     molecular weight can be calculated at each level

      CALL PROMPT('Enter existing filename : ')
      READ(*,10)IPFILE
10    FORMAT(A)
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')


      CALL PROMPT('Enter output filename : ') 
      READ(*,10)OPFILE
      OPEN(UNIT=2,FILE=OPFILE,STATUS='UNKNOWN')

C     First skip header
54    READ(1,1)BUFFER
1     FORMAT(A)
      IF(BUFFER(1:1).EQ.'#') THEN
       WRITE(2,1)BUFFER
       GOTO 54
      ENDIF
      READ(BUFFER,*)AMFORM
      WRITE(2,*)AMFORM
      IF(AMFORM.EQ.0)THEN
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
        WRITE(2,501)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
      ELSE
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
        WRITE(2,511)IPLANET,LATITUDE,NPRO,NVMR
      ENDIF
501   FORMAT(1X,I3,1X,F6.2,1X,I3,I3,F8.3)
511   FORMAT(1X,I3,1X,F6.2,1X,I3,I3)

      DO 20 I=1,NVMR
         READ(1,*)ID(I),ISO(I)
         WRITE(2,502)ID(I),ISO(I)
20    CONTINUE
502   FORMAT(1X,I3,I5)

C     reading the first block of profiles
      READ(1,*)
      N=MIN(NVMR,3)
C     N is the maximum VMR which can be read in from the next block
      DO 30 I=1,NPRO
         READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,N)
30    CONTINUE
C     reading in additional blocks if any
C     N VMR profiles have been read in so far
33    IF(NVMR.GT.N)THEN
       READ(1,*)
C      profiles up to VMR(?,K) to be read from this block
       K=MIN(NVMR,(N+6))
       DO 32 I=1,NPRO
        READ(1,*)(VMR(I,J),J=N+1,K)
32     CONTINUE
       N=K
       GOTO 33
      END IF
      CLOSE(UNIT=1)

      WRITE(2,504)(I,I=1,NVMR)
504   FORMAT(1X,' height (km) ',' press (atm) ','  temp (K)   ',
     1  40(' VMR gas',I3,2X))
      DO 505 I=1,NPRO
        WRITE(2,506)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
505   CONTINUE
506   FORMAT(1X,F13.3,E13.5,F13.4,40(E13.5))

      CLOSE(UNIT=2)

      END


