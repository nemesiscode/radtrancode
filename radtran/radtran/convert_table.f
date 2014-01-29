      PROGRAM CONVERT_TABLE
C     ********************************************************************
C     Program to convert Radtran .kta binary k-coefficient files into more
C     easily read ASCII .par format
C
C     Pat Irwin      31/1/96
C     ********************************************************************

C     note dbcom defines the linedata base variables. it is not normally stored
C     in the same directory as the rest of the code
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/bincom.f'
C-----------------------------------------------------------------------------
      INTEGER LUN,LOOP,LUN0,LUN1,NW
      PARAMETER (LUN=2,LUN0=30,LUN1=36)
C     MAXOUT the maximum number of output points
      INTEGER IREC,IREC0,I,ITAB
      CHARACTER*100 KTAFIL,PARFIL
      REAL TE1,TEMP1(MAXK),PRESS1(MAXK)
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG)


      WRITE(*,*)'Enter number of 4-byte words per record.'
      CALL PROMPT('old files = 2, new files = 1 : ')  
      READ*,NW
      IF(NW.LT.1.OR.NW.GT.2)THEN
       PRINT*,'NW must be 1 or 2 only!!!'
       STOP
      ENDIF

      IRECL = NW*ISYS()
      print*,'irecl = ',irecl

      CALL PROMPT('Enter input ktable filename : ')
      READ(5,23)OPFILE
23    FORMAT(A)
      CALL FILE(OPFILE,KTAFIL,'kta')

      CALL FILE(OPFILE,PARFIL,'par')

      OPEN(LUN1,FILE=PARFIL,STATUS='UNKNOWN')

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1RECL=IRECL)
      READ(LUN0,REC=1)IREC0
      READ(LUN0,REC=2)NPOINT
      READ(LUN0,REC=3)VMIN
      READ(LUN0,REC=4)DELV
      READ(LUN0,REC=5)FWHM
      READ(LUN0,REC=6)NP
      READ(LUN0,REC=7)NT
      READ(LUN0,REC=8)NG
      READ(LUN0,REC=9)IDGAS(1)
      READ(LUN0,REC=10)ISOGAS(1)
      IREC=11
      NGAS=1

      IF(NG.NE.10.AND.NG.NE.20)THEN
       PRINT*,'Convert_table is currently limited to NG=10 or NG=20'
       PRINT*,'Aborting'
       STOP
      ENDIF
      WRITE(LUN1,*)'VMIN,DELV,FWHM = ',VMIN,DELV,FWHM
      print*,'VMIN,DELV,FWHM = ',VMIN,DELV,FWHM
      WRITE(LUN1,*)'NPOINT,VMAX = ',NPOINT,VMIN + (NPOINT-1)*DELV
      print*,'NPOINT,VMAX = ',NPOINT,VMIN + (NPOINT-1)*DELV
      WRITE(LUN1,*)'Gas ID,ISO : ',IDGAS(1),ISOGAS(1)
      print*,'NP,NT,NG = ',NP,NT,NG
      print*,'ID, ISO = ',IDGAS(1),ISOGAS(1)
      print*,' '
      DO 299 J=1,NG
       READ(LUN0,REC=IREC)G_ORD(J)
       IREC=IREC+1
299   CONTINUE
      print*,'G - ordinates, weights'
      WRITE(LUN1,*)'Number of g-ordinates: ',NG
      WRITE(LUN1,*)'G_ORD, DEL_G'
      DO 399 J=1,NG
       READ(LUN0,REC=IREC)DEL_G(J)
       print*,g_ord(j),del_g(j)
       WRITE(LUN1,*)G_ORD(J),DEL_G(J)
       IREC=IREC+1
399   CONTINUE
      IREC=11 + 2*NG + 2
      print*,'Pressures : '
      WRITE(LUN1,*)'Number of pressures : ',NP
      WRITE(LUN1,*)'Pressures (atm)'
      DO 301 J=1,NP
       READ(LUN0,REC=IREC)PRESS1(J)
       print*,press1(j)
       WRITE(LUN1,*)press1(j)
       IREC=IREC+1
301   CONTINUE
      print*,'Temperatures : '
      WRITE(LUN1,*)'Number of temperatures : ',NT
      WRITE(LUN1,*)'Temperatures (K)'
      DO 302 J=1,NT
       READ(LUN0,REC=IREC)TEMP1(J)
       print*,temp1(j)
       WRITE(LUN1,*)temp1(j)
       IREC=IREC+1
302   CONTINUE

      IREC=IREC0

      CALL PROMPT('Data tab. by wavenumber (0) or wavelength (1)? : ')
      READ*,ITAB

      WRITE(LUN1,*)'Absorption coefficients  (km.amagat)-1 : '

      DO 1000 I=1,NPOINT

      VV = VMIN + (I-1)*DELV


      IREC=IREC0+NP*NT*NG*(I-1)

       DO 30 K=1,NT

         TE1=TEMP1(K)
         IF(ITAB.EQ.0)THEN
          WRITE(LUN1,996)VV,TE1
         ELSE
          WRITE(LUN1,997)VV,TE1
         ENDIF 
996      FORMAT('          RESULTS FOR',F7.1,' CM-1,',F6.1,
     &' DEGREES K')
997      FORMAT('          RESULTS FOR',F7.4,' MICRON,',F6.1,
     &' DEGREES K')
        WRITE(LUN1,*)' '

        WRITE(LUN1,*)'Pressure, and fitted K-Coefficients'

        DO 20 J=1,NP
         PE1 = PRESS1(J)
         DO 40 LOOP=1,NG
          READ(LUN0,REC=IREC)K_G(LOOP)
          IREC=IREC+1
40       CONTINUE
         WRITE(LUN1,707) PE1,(28650.0*K_G(LOOP),LOOP=1,5)
         WRITE(LUN1,708) (28650.0*K_G(LOOP),LOOP=6,10)
         IF(NG.EQ.20)THEN
          WRITE(LUN1,708) (28650.0*K_G(LOOP),LOOP=11,15)
          WRITE(LUN1,708) (28650.0*K_G(LOOP),LOOP=16,20)
         ENDIF
20      CONTINUE
30     CONTINUE

1000  CONTINUE

707    FORMAT(E11.4,5(E13.6))
708    FORMAT(11X,5(E13.6))



      END



