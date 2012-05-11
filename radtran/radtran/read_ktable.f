      SUBROUTINE READ_KTABLE(KTAFIL,LUN0,NPOINT,VMIN,DELV,FWHM,
     1IDGAS,ISOGAS,PRESS1,TEMP1,NP,NT,G_ORD,DEL_G,NG,IREC0)
C     $Id: read_ktable.f,v 1.4 2011-06-17 15:40:27 irwin Exp $
C---------------------------------------------------------------------------
C_TITLE:  READ_KTABLE
C
C_ARGS:
C
C_KEYS:
C
C_DESCR:  reads an absorption coefficient look-up table
C
C_FILES:  UNIT LUN: output file
C
C_CALLS:
C
C_BUGS:   
C
C_HIST:   10feb95 PGJI ORIGINAL VERSION
C---------------------------------------------------------------------------
      INCLUDE '../includes/arrdef.f'
      INTEGER PINDEX,CP,CT,LOOP,LUN0
C     MAXOUT the maximum number of output points
      INTEGER NOUT,LIMLAY,INCDIM,MAXPT,ICON,IREC,IREC0,NSKIP,I
      real t1,t2,tmp(2)
      CHARACTER*100 DRVFIL,KTAFIL
      LOGICAL BIT
      REAL PP1(MAXK),P1,TE1,TEMP1(MAXK),PRESS1(MAXK)
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),TABLE(MAXK,MAXK,MAXG)
      INTEGER IDGAS,ISOGAS,IDGAS1,ISOGAS1

      REAL VMIN1,DELV1,FWHM1
      INTEGER NPOINT1


      IRECL=ISYS()

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1RECL=IRECL)
      READ(LUN0,REC=1)IREC0
      READ(LUN0,REC=2)NPOINT1
      READ(LUN0,REC=3)VMIN1
      READ(LUN0,REC=4)DELV1
      READ(LUN0,REC=5)FWHM1
      READ(LUN0,REC=6)NP
      READ(LUN0,REC=7)NT
      READ(LUN0,REC=8)NG
      READ(LUN0,REC=9)IDGAS1
      READ(LUN0,REC=10)ISOGAS1

      PRINT*,'IREC0 = ',IREC0
      IF(IDGAS1.NE.IDGAS)THEN
       PRINT*,'Read_ktable. IDGAS1 <> IDGAS'
       PRINT*,'Table has : ',IDGAS1,ISOGAS1
       PRINT*,'Expecting : ',IDGAS,ISOGAS
       STOP
      END IF

      IF(ISOGAS1.NE.ISOGAS)THEN
       PRINT*,'Read_ktable. ISOGAS1 <> ISOGAS'
       STOP
      END IF

      IF(DELV1.NE.DELV)THEN
       PRINT*,'Read_ktable. DELV1 <> DELV'
       STOP
      END IF

      IF(FWHM1.NE.FWHM)THEN
       PRINT*,'Read_ktable. FWHM1 <> FWHM'
       STOP
      END IF

      VMAX = VMIN + (NPOINT-1)*DELV
      VMAX1 = VMIN1 + (NPOINT1-1)*DELV1

      IF(VMIN.LT.VMIN1)THEN
       PRINT*,'Read_ktable. VMIN < VMIN1'
       STOP
      END IF

      IF(VMAX.GT.VMAX1)THEN
       PRINT*,'Read_ktable. VMAX > VMAX1'
       STOP
      END IF

      IREC=11
      DO 299 J=1,NG
       READ(LUN0,REC=IREC)G_ORD(J)
       IREC=IREC+1
299   CONTINUE
      DO 399 J=1,NG
       READ(LUN0,REC=IREC)DEL_G(J)
       IREC=IREC+1
399   CONTINUE
      IREC=11 + 2*NG + 2
      DO 301 J=1,NP
       READ(LUN0,REC=IREC)PRESS1(J)
       PRESS1(J)=LOG(PRESS1(J))
       IREC=IREC+1
301   CONTINUE
      DO 302 J=1,NT
       READ(LUN0,REC=IREC)TEMP1(J)
       IREC=IREC+1
302   CONTINUE
 
      N1 = 1 + INT((VMIN-VMIN1)/DELV)

      IREC0=IREC0+NP*NT*NG*(N1-1)
      PRINT*,VMIN,VMIN1,N1,IREC0

      RETURN
      END


      SUBROUTINE INTERP_KTABLE(PRESS1,TEMP1,NP,NT,NG,LUN0,IREC0,
     1PIN,TIN,V,VMIN,DELV,K_G)
C     *******************************************************************
C     *******************************************************************
      INCLUDE '../includes/arrdef.f'
      REAL TEMP1(MAXK),PRESS1(MAXK),PIN,TIN
      INTEGER NP,NT,NG,LUN0,IREC0
      REAL V,VMIN,DELV,K_G(MAXG)

      INTEGER N1,IREC,CT,CP,J,K,LOOP,PINDEX
      REAL TABLE(20,20,20),P1,T1
      REAL Y1,Y2,Y3,Y4,U,T

C     ------------------------------------------------------------------

      N1 = 1 + INT((V - VMIN)/DELV)

      IREC=IREC0+NP*NT*NG*(N1-1)

       DO 20 J=1,NP
        DO 30 K=1,NT

         DO 40 LOOP=1,NG
          READ(LUN0,REC=IREC)TABLE(J,K,LOOP)
          IREC=IREC+1
40       CONTINUE

30      CONTINUE
20     CONTINUE
10    CONTINUE


      PMAX=PRESS1(NP)
      PMIN=PRESS1(1)
      P1=LOG(PIN)
      T1=TIN
      IF(P1.LT.PMIN)THEN
C       Print*,'Warning P<Pmin. Setting to Pmin',PMIN
       P1=PMIN
      END IF
      IF(P1.GT.PMAX)THEN
C       Print*,'Warning P>Pmax. Setting to Pmax',PMAX
       P1=PMAX
      END IF
      TMAX=TEMP1(NT)
      TMIN=TEMP1(1)
      IF(TIN.LT.TMIN)THEN
C       Print*,'Warning T<Tmin. Setting to Tmin',TMIN
       T1=TMIN
      END IF
      IF(TIN.GT.TMAX)THEN
C       Print*,'Warning T>Tmax. Setting to Tmax',TMAX
       T1=TMAX
      END IF

      CP = PINDEX(PRESS1,NP,P1)

      CT = PINDEX(TEMP1,NT,T1)

      DO 80 LOOP=1,NG
       Y1=TABLE(CP,CT,LOOP)
       Y2=TABLE(CP+1,CT,LOOP)
       Y3=TABLE(CP+1,CT+1,LOOP)
       Y4=TABLE(CP,CT+1,LOOP)
       T=(P1-PRESS1(CP))/(PRESS1(CP+1)-PRESS1(CP))
       U=(T1-TEMP1(CT))/(TEMP1(CT+1)-TEMP1(CT))
       K_G(LOOP)=(1.0-T)*(1.0-U)*Y1 + T*(1.0-U)*Y2 + T*U*Y3 + 
     1  (1.0-T)*U*Y4

80    CONTINUE

      RETURN

      END



      INTEGER FUNCTION PINDEX(X,NX,X1)
      INTEGER NX
      REAL X(NX),X1
      IF(X1.GT.X(NX))THEN
C       PRINT*,'Warning in INDEX. X1>XMAX',X1,X(NX)
       X1=X(NX)
      END IF
      IF(X1.LT.X(1))THEN
C       PRINT*,'Warning in INDEX. X1>XMIN',X1,X(1)
       X1=X(1)
      END IF

      J=0

      DO 10 I=1,NX-1
       IF(X1.GE.X(I))THEN
        J=I
       ELSE
        GOTO 20
       END IF
10    CONTINUE
20    CONTINUE

      PINDEX=J
      RETURN
      END
