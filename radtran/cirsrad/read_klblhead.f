************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE READ_KLBLHEAD
C
C	Opens LTA file, reads in header,and returns record number for 
C	the beginning of the range of interest as well as temp and 
C	pressure ranges used. This differs from READ_KTABLE in that the
C	latter checks that the input file matches gas and isotope IDs
C	deemed required by the calling routine.
C
C	PGJI 		10/02/95
C	A.L.Weir	03/03/96
C	PGJI		15/11/16 modified from read_khead.f
C
C-----------------------------------------------------------------------

      SUBROUTINE READ_KLBLHEAD(KTAFIL,LUN0,NPOINT,VMIN,DELV,
     1		IDGAS,ISOGAS,PRESS,TEMP,TEMP2,NP,NT,IREC0)

	implicit none
        include '../includes/arrdef.f'
	integer		LUN0, NPOINT, IDGAS, ISOGAS, NP, NT, IREC0,
     1			IREC, N1, IRECL, J, ISYS, I
	real		VMIN, DELV, PRESS(MAXK), TEMP(MAXK), 
     1                  TEMP2(MAXK,MAXK)
	character*100	KTAFIL,IPFILE

      IRECL=ISYS()

      WRITE(*,41)KTAFIL
41    FORMAT('Read_klblhead. Opening : ',A) 

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1RECL=IRECL)

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
      READ(LUN0,REC=5)NP
      READ(LUN0,REC=6)NT
      READ(LUN0,REC=7)IDGAS
      READ(LUN0,REC=8)ISOGAS


C-----------------------------------------------------------------------
C
C	Read in g ordinates, deltag, pressures and temperatures.
C
C-----------------------------------------------------------------------

      IREC=9
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
 
      RETURN
      END

************************************************************************
************************************************************************
