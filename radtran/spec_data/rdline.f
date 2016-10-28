      SUBROUTINE RDLINE(BUFFER)
C     $Id: rdline.f,v 1.9 2011-06-17 14:54:42 irwin Exp $
C***********************************************************************
C_TITL:	RDLINE.f
C
C_DESC:	Reads line data parameters from a text buffer. Uses an internal
C	READ to read line data parameters from a text buffer. The format
C	of the READ is determined by DBFORM. Any variables not included in
C	the data base format are set to null or safe values.
C
C_ARGS:	Input variable:
C	BUFFER	CHARA*(*)	Text buffer containing the line data record.
C				CHARACTER*(*) declares an incoming
C				character variable whose length is   
C				unknown.
C
C_FILE:	No files openned.
C
C_CALL: No calls made.
C
C_HIST:	4feb91	SBC	Original version
C***************************** VARIABLES *******************************

      CHARACTER*(*) BUFFER

C ../includes/dbcom.f stores the linedata base variables (e.g. RELABU).
      INCLUDE '../includes/dbcom.f' 

      CHARACTER*3 ATNW
      INTEGER EXP,EXP1,IPOW,IPOW1
      DOUBLE PRECISION CC
      PARAMETER (CC=10.)
      COMMON /OVERFLOW/IPOW,IPOW1

C******************************** CODE *********************************

cc      WRITE(*,*)' RDLINE.f :: DBFORM, DBRECL = ',DBFORM,DBRECL
      DBFORM=0
      IF(DBFORM.EQ.0)THEN
C=======================================================================
C
C     HITRAN (either of 100-, 112-character or 160-character formats)
C
C=======================================================================
        IF(DBRECL.EQ.100)THEN
          READ(BUFFER,100,ERR=99)LNID,LNISO,LNWAVE,LNSTR,EXP,LNPROB,
     1    EXP1,LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,LNUGQI,LNLGQI,LNULQ,
     2    LNLLQ,LNACC,LNREF
100       FORMAT(I2,I1,F12.6,F6.3,1X,I3,F6.3,1X,I3,F5.4,F5.4,F10.4,
     1    F4.2,F8.6,I3,I3,A9,A9,3I1,3I2)
          LDOUBV= 0.0
        ELSE IF(DBRECL.EQ.112)THEN
          READ(BUFFER,101,ERR=99)LNID,LNISO,LNWAVE,LNSTR,EXP,LNPROB,
     1    EXP1,LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,LNUGQI,LNLGQI,LNULQ,
     2    LNLLQ,LNACC,LNREF,LDOUBV
101       FORMAT(I2,I1,F12.6,F6.3,1X,I3,F6.3,1X,I3,F5.4,F5.4,F10.4,
     1    F4.2,F8.6,I3,I3,A9,A9,3I1,3I2,F12.7)
        ELSE IF(DBRECL.EQ.160)THEN
C         HITEMP160 is not quite the same as HITRAN160. Need to add a 
C         catch to spot the difference.
          IF(BUFFER(21:21).EQ.'E')THEN
C          HITEMP160
           READ(BUFFER,102,ERR=99)LNID,LNISO,LNWAVE,LNSTR,EXP,LNEINA,
     1     EXP1,LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,LNUGQI04,LNLGQI04,
     2     LNULQ04,LNLLQ04,LNACC04,LNREF04,LNFLAG,UWGHT,LWGHT
          ELSE
C          HITRAN160
           READ(BUFFER,104,ERR=99)LNID,LNISO,LNWAVE,LNSTR,EXP,LNEINA,
     1     EXP1,LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,LNUGQI04,LNLGQI04,
     2     LNULQ04,LNLLQ04,LNACC04,LNREF04,LNFLAG,UWGHT,LWGHT
          ENDIF
102       FORMAT(I2,I1,F12.6,F5.3,1X,I4,F6.3,1X,I3,F5.4,F5.4,F10.4,
     1    F4.2,F8.6,A15,A15,A15,A15,6I1,6I2,A1,F7.1,F7.1)
104       FORMAT(I2,I1,F12.6,F6.3,1X,I3,F6.3,1X,I3,F5.4,F5.4,F10.4,
     1    F4.2,F8.6,A15,A15,A15,A15,6I1,6I2,A1,F7.1,F7.1)
        ELSE IF(DBRECL.EQ.52)THEN
          IF(BUFFER(21:21).EQ.'E')THEN
          READ(BUFFER,107,ERR=99)LNID,LNISO,LNWAVE,LNSTR,EXP,
     1    LNLSE,LNWIDA,LNTDEP,LNWIDS,LNTDEPS
          ELSE
          READ(BUFFER,106,ERR=99)LNID,LNISO,LNWAVE,LNSTR,EXP,
     1    LNLSE,LNWIDA,LNTDEP,LNWIDS,LNTDEPS
          ENDIF
106       FORMAT(I2,I1,F12.6,F6.3,1X,I3,F10.4,F5.4,F3.2,F6.4,F3.2)
107       FORMAT(I2,I1,F12.6,F5.3,1X,I4,F10.4,F5.4,F3.2,F6.4,F3.2)
        ELSE
          WRITE(*,*)' RDLINE.f :: HITRAN format not recognised.'
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' Accepted DBRECL for DBFORM = 0, either 100,112.'
          WRITE(*,*)' DBFORM, DBRECL = ',DBFORM,DBRECL
          STOP
        ENDIF
        
       IF(DBRECL.EQ.52)THEN
        LNEINA=0.0
        LNPROB= 0.0
        EXP1= 1
        LNPSH= 0.0
        LNUGQI04= ''
        LNLGQI04= ''
        LNULQ04= ' '
        LNLLQ04= ' '
        LNACC04=0
        LNFLAG=''
        UWGHT=0
        LWGHT=0
       ENDIF
      ELSE IF(DBFORM.EQ.1)THEN
C=======================================================================
C
C     GEISA (either of 80-,120- or 211-character formats)
C
C=======================================================================
        IF(DBRECL.EQ.82)THEN
          READ(BUFFER,210,ERR=99)LNWAVE,LNSTR,EXP,LNWIDA,LNLSE,LNTDEP,
     1    LNISO,LNID
210       FORMAT(F10.3,F6.3,1X,I3,F5.3,F10.3,36X,F4.2,I4,I3)
C NOTE: above assumes that GEISA version with temp dependence is
C used. Also ignores quantum numbers in this format
          LNPROB= 0.0
          EXP1= 1
          LNWIDS= 0.0
          LNPSH= 0.0
          LNUGQI= 0
          LNLGQI= 0
          LNULQ= ' '
          LNLLQ= ' '
          DO I=1,3
            LNACC(I)= 0
            LNREF(I)= 0
          ENDDO
          LNTDEPS= LNTDEP
          LDOUBV= 0.0
        ELSE IF(DBRECL.EQ.80)THEN
          READ(BUFFER,200,ERR=99)LNWAVE,LNSTR,EXP,LNWIDA,LNLSE,ATNW,
     1    LNISO,LNID
200       FORMAT(F10.3,F6.3,1X,I3,F5.3,F10.3,35X,A3,I4,I3)
C NOTE: above assumes that GEISA version with temp dependence is
C used. Also ignores quantum numbers in this format
          LNPROB= 0.0
          EXP1= 1
          LNWIDS= 0.0
          LNPSH= 0.0
          LNUGQI= 0
          LNLGQI= 0
          LNULQ= ' '
          LNLLQ= ' '
          DO I=1,3
            LNACC(I)= 0
            LNREF(I)= 0
          ENDDO
          LNTDEP= 0.0
          LNTDEPS= 0.0
          LDOUBV= 0.0
          IF(ATNW(1:1).EQ.'0'.OR.ATNW(1:1).EQ.'.')THEN
            READ(ATNW,*)LNTDEP
            LNTDEPS= LNTDEP
          ENDIF
C Hack to cure line strength bug with CH3D in Geisa '84
          IF (LNID.EQ.23) THEN
            LNSTR= LNSTR*RELABU(3,6)
          ENDIF


        ELSE IF(DBRECL.EQ.120)THEN
          READ(BUFFER,205,ERR=99)LNWAVE,LNSTR,EXP,LNWIDA,LNLSE,LNTDEP,
     1    LNISO,LNID,LNPROB,EXP1,LNWIDS,LNPSH,LNACC,LNREF
205       FORMAT(F10.3,F6.3,1X,I3,F5.3,F10.3,36X,F4.2,I4,I3,6X,
     1    F6.3,1X,I3,F5.4,F8.6,3I1,3I2)
          LNTDEPS= LNTDEP
          LNUGQI= 0
          LNLGQI= 0
          LNULQ= ' '
          LNLLQ= ' '
          LDOUBV= 0.0
        ELSE IF(DBRECL.EQ.211)THEN
          READ(BUFFER,206,ERR=99)LNWAVE,LNSTR,EXP,LNWIDA,LNLSE,LNTDEP,
     1    LNISO,LNID,LNPROB,EXP1,LNWIDS,LNPSH,LNACC,LNREF
206       FORMAT(F12.6,F7.4,1X,I3,F6.4,F10.4,36X,F4.2,I3,I3,6X,
     1    F6.3,1X,I3,F5.4,F8.6,3I1,3I2)
          LNTDEPS= LNTDEP
C          PRINT*,LNWAVE,LNSTR,EXP,LNWIDA,LNLSE,LNTDEP,
C     1    LNISO,LNID,LNPROB,EXP1,LNWIDS,LNPSH,LNACC,LNREF
          LNUGQI= 0
          LNLGQI= 0
          LNULQ= ' '
          LNLLQ= ' '
          LDOUBV= 0.0 
        ELSE
    	  WRITE(*,*)' RDLINE.f :: GEISA format not recognised.'
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' Accepted DBRECL for DBFORM = 1, either 80,82,120.'
          WRITE(*,*)' DBFORM, DBRECL = ',DBFORM,DBRECL
          STOP
        ENDIF
      ELSE IF(DBFORM.EQ.2)THEN
C HITRAN with temperature coefficient of self-broadening in the transition
C probability column.
         LNPROB= 0.0
         READ(BUFFER,300,ERR=99)LNID,LNISO,LNWAVE,LNSTR,EXP,LNTDEPS,
     1   LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,
     2   LNUGQI,LNLGQI,LNULQ,LNLLQ,LNACC,LNREF
300      FORMAT(I2,I1,F12.6,F6.3,1X,I3,E10.3,F5.4,F5.4,F10.4,
     1   F4.2,F8.6,I3,I3,A9,A9,3I1,3I2)
      ELSE
        WRITE(*,*)' RDLINE.f :: Invalid data format. Stopping program.'
        STOP
      ENDIF

C Scaling strengths by 1.E47 to avoid underflow and including exponent
C Old versions used Avagadros constant as the scaling factor but new data
C bases include many very weak lines
       if(EXP.LE.(-14))then
        LNSTR= LNSTR*CC**(EXP+47)
       elseif(EXP.LT.(-14).AND.EXP.GE.(-20))THEN
        print*, 'Note: LNSTR EXP = ', EXP
        print*, 'Possibly suspect!'
        LNSTR= LNSTR*CC**(EXP+47)
       else
        LNSTR=0
        print*, 'LNSTR exp > -14, setting to zero: rdline.f'
        print*, 'EXP = ', EXP, LNWAVE
       endif
C      print*,LNSTR,EXP
      
C      print*,'A',LNSTR
      IF(DBRECL.EQ.160)THEN
       LNEINA= LNEINA*CC**EXP1
      ELSE
       LNPROB= LNPROB*CC**EXP1
      ENDIF

      IF(LNLSE.LT.0.0)THEN
C       WRITE(*,*)' '
C       WRITE(*,*)' RDLINE.f :: ELIN < 0 '
C       WRITE(*,*)' No temperature variation of line strength'
C       WRITE(*,103)BUFFER
C       DO 13 I=1,DBRECL,20
C        WRITE(*,*)I,'!',BUFFER(I:I+19),'@'
C13     CONTINUE
C       WRITE(*,*)' Setting ELIN to zero'
       LNLSE=0.0
      ENDIF

      RETURN

99    CONTINUE

C=======================================================================
C
C	For debugging purposes ...
C
C=======================================================================
      WRITE(*,*)' '
      WRITE(*,*)' RDLINE.f :: Could not read buffer.'
      WRITE(*,*)' BUFFER = '
      WRITE(*,103)BUFFER
103   FORMAT(A)
C      DO 13 I=1,DBRECL,20
C        WRITE(*,*)I,'!',BUFFER(I:I+19),'@'
C13    CONTINUE

      print*,LNID,LNISO,LNWAVE,LNSTR,EXP
      print*,LNEINA,EXP1,LNWIDA,LNWIDS
      print*,LNLSE,LNTDEP,LNPSH

C End of debugging lines

      LNID= 0
C Setting LNWAVE to zero so that any search routines finding leading 
C comments behave correctly
      LNWAVE= 0.0

      RETURN

      END
