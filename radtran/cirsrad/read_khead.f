************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE READ_KHEAD
C
C	Opens KTA file, reads in header,and returns record number for 
C	the beginning of the range of interest as well as temp and 
C	pressure ranges used. This differs from READ_KTABLE in that the
C	latter checks that the input file matches gas and isotope IDs
C	deemed required by the calling routine.
C
C	PGJI 		10/02/95
C	A.L.Weir	03/03/96
C
C-----------------------------------------------------------------------

      SUBROUTINE READ_KHEAD(KTAFIL,LUN0,NPOINT,VMIN,DELV,FWHM,VCEN,
     1		IDGAS,ISOGAS,PRESS,TEMP,TEMP2,NP,NT,G_ORD,DEL_G,NG,
     2          IREC0)

	implicit none
        include '../includes/arrdef.f'
	integer		LUN0, NPOINT, IDGAS, ISOGAS, NP, NT, NG,
     1			IREC, N1, IRECL, J, ISYS, I, IREC0
	real		VMIN, DELV, FWHM, PRESS(MAXK), TEMP(MAXK), 
     1			G_ORD(MAXG), DEL_G(MAXG), VCEN(MAXBIN),
     2                  TEMP2(MAXK,MAXK)
	character*200	KTAFIL
        integer idiag,iquiet
        common/diagnostic/idiag,iquiet
   
      IRECL=ISYS()

C      WRITE(*,41)KTAFIL
41    FORMAT('Read_khead. Opening : ',A) 

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1RECL=IRECL, ACTION='READ')

      READ(LUN0,REC=1)IREC0
      READ(LUN0,REC=2)NPOINT
      IF(NPOINT.GT.MAXBIN)THEN
       PRINT*,'Problem in read_khead.f. Number of points in k-table'
       PRINT*,'is greater than the dimension, MAXBIN,of the VCEN array'
       PRINT*,NPOINT,MAXBIN
       STOP
      ENDIF
      READ(LUN0,REC=3)VMIN
      READ(LUN0,REC=4)DELV
      READ(LUN0,REC=5)FWHM
      READ(LUN0,REC=6)NP
      READ(LUN0,REC=7)NT
      READ(LUN0,REC=8)NG
      READ(LUN0,REC=9)IDGAS
      READ(LUN0,REC=10)ISOGAS

C-----------------------------------------------------------------------
C
C	Read in g ordinates, deltag, pressures and temperatures.
C
C-----------------------------------------------------------------------

C     Check to see if NG is within limits
      IF(NG.GT.MAXG-1)THEN
       print*,'Error in read_khead.f: NG is greater than allowed'
       print*,'maximum value MAXG, defined in arrdef.f'
       print*,NG,MAXG
       stop
      ENDIF
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
       READ(LUN0,REC=IREC)PRESS(J)
       PRESS(J)=LOG(PRESS(J))
       IREC=IREC+1
301   CONTINUE
      N1=ABS(NT)
      IF(NT.LT.0)THEN
       DO 307 I=1,NP
        DO 308 J=1,N1
        READ(LUN0,REC=IREC)TEMP2(I,J)
        IREC=IREC+1 
308     CONTINUE
307    CONTINUE
      ELSE
       DO 302 J=1,NT
        READ(LUN0,REC=IREC)TEMP(J)
        IREC=IREC+1
302    CONTINUE
      ENDIF
      IF(DELV.LE.0.0)THEN
       DO 303 J=1,NPOINT   
       READ(LUN0,REC=IREC)VCEN(J)
       IREC = IREC + 1
303    CONTINUE
      ELSE
       DO 304 J=1,NPOINT  
        VCEN(J)=VMIN+FLOAT(J-1)*DELV
304    CONTINUE
      ENDIF

 
C-----------------------------------------------------------------------
C
C	Finally, return the record number for the beginning of the
C	wavenumber region of interest. 
C
C-----------------------------------------------------------------------

      RETURN
      END

************************************************************************
************************************************************************
