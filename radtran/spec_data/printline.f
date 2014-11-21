      SUBROUTINE PRINTLINE(BUFFER,NFORM,NEWRECL)
C***********************************************************************
C_TITL:	PRINTLINE.f
C
C_DESC:	prints out line data parameters to buffer in various GEISA and HITRAN
C       formats. 
C	The format of the output is determined by NFORM and NEWRECL. 
C
C_ARGS:	Input variables:
C	NFORM	INTEGER	0=HITRAN,1=GEISA,2=OTHER
C	NEWRECL	INTEGER	Record length
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
C	11jun12	PGJI	Modified from RDLINE to print out rather than read in.
C
C***************************** VARIABLES *******************************

      CHARACTER*(*) BUFFER

C ../includes/dbcom.f stores the linedata base variables (e.g. RELABU).
      INCLUDE '../includes/dbcom.f' 

      INTEGER NEWRECL,NFORM
      CHARACTER*1 ANS
      REAL L1

C******************************** CODE *********************************

C     Scale lines.
      L1 = LNSTR*1e-20
      L1 = L1*1e-27


      IF(NFORM.EQ.0)THEN
C=======================================================================
C
C     HITRAN (either of 100-, 112-character or 160-character formats)
C
C=======================================================================
        IF(NEWRECL.EQ.100)THEN
          WRITE(BUFFER,100,ERR=99)LNID,LNISO,LNWAVE,L1,LNPROB,
     1    LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,LNUGQI,LNLGQI,LNULQ,
     2    LNLLQ,LNACC,LNREF
100       FORMAT(I2,I1,F12.6,E10.3,E10.3,F5.4,F5.4,F10.4,
     1    F4.2,F8.6,I3,I3,A9,A9,3I1,3I2)
        ELSE IF(NEWRECL.EQ.112)THEN
          WRITE(BUFFER,101,ERR=99)LNID,LNISO,LNWAVE,L1,LNPROB,
     1    LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,LNUGQI,LNLGQI,LNULQ,
     2    LNLLQ,LNACC,LNREF,LDOUBV
101       FORMAT(I2,I1,F12.6,E10.3,E10.3,F5.4,F5.4,F10.4,
     1    F4.2,F8.6,I3,I3,A9,A9,3I1,3I2,F12.7)
        ELSE IF(NEWRECL.EQ.160)THEN
          WRITE(BUFFER,102,ERR=99)LNID,LNISO,LNWAVE,L1,LNEINA,
     1    LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,LNUGQI04,LNLGQI04,
     2    LNULQ04,LNLLQ04,LNACC04,LNREF04,LNFLAG,UWGHT,LWGHT
102       FORMAT(I2,I1,F12.6,E10.3,E10.3,F5.4,F5.4,F10.4,
     1    F4.2,F8.6,A15,A15,A15,A15,6I1,6I2,A1,F7.1,F7.1)
        ELSE
          WRITE(*,*)' PRINTLINE.f :: HITRAN format not recognised.'
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' Accepted DBRECL for NFORM = 0, either 100,112.'
          WRITE(*,*)' NFORM, DBRECL = ',NFORM,DBRECL
          STOP
        ENDIF
 
      ELSE IF(NFORM.EQ.1)THEN
C=======================================================================
C
C     GEISA (either of 80-,120- or 211-character formats)
C
C=======================================================================
        IF(NEWRECL.EQ.82)THEN
          WRITE(BUFFER,210,ERR=99)LNWAVE,L1,LNWIDA,LNLSE,LNTDEP,
     1    LNISO,LNID
210       FORMAT(F10.3,E10.3,F5.3,F10.3,36X,F4.2,I4,I3)


C NOTE: above assumes that GEISA version with temp dependence is
C used. Also ignores quantum numbers in this format

        ELSE IF(NEWRECL.EQ.80)THEN

          IF (LNID.EQ.23) THEN
            L1= L1/RELABU(3,6)
          ENDIF

          WRITE(BUFFER,200,ERR=99)LNWAVE,L1,EXP,LNWIDA,LNLSE,LNTDEP,
     1    LNISO,LNID
200       FORMAT(F10.3,E10.3,F5.3,F10.3,35X,F3.2,I4,I3)


        ELSE IF(NEWRECL.EQ.120)THEN
          WRITE(BUFFER,205,ERR=99)LNWAVE,L1,LNWIDA,LNLSE,LNTDEP,
     1    LNISO,LNID,LNPROB,LNWIDS,LNPSH,LNACC,LNREF
205       FORMAT(F10.3,E10.3,F5.3,F10.3,36X,F4.2,I4,I3,6X,
     1    E10.3,F5.4,F8.6,3I1,3I2)
        ELSE IF(NEWRECL.EQ.211)THEN
          WRITE(BUFFER,206,ERR=99)LNWAVE,L1,LNWIDA,LNLSE,LNTDEP,
     1    LNISO,LNID,LNPROB,LNWIDS,LNPSH,LNACC,LNREF
206       FORMAT(F12.6,F11.4,F6.4,F10.4,36X,F4.2,I3,I3,6X,
     1    E10.3,F5.4,F8.6,3I1,3I2)
        ELSE
    	  WRITE(*,*)' PRINTLINE.f :: GEISA format not recognised.'
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' Accepted DBRECL for NFORM = 1, either 80,82,120.'
          WRITE(*,*)' NFORM, DBRECL = ',NFORM,DBRECL
          STOP
        ENDIF
      ELSE IF(NFORM.EQ.2)THEN
C HITRAN with temperature coefficient of self-broadening in the transition
C probability column.
         LNPROB= 0.0
         WRITE(BUFFER,300,ERR=99)LNID,LNISO,LNWAVE,L1,LNTDEPS,
     1   LNWIDA,LNWIDS,LNLSE,LNTDEP,LNPSH,
     2   LNUGQI,LNLGQI,LNULQ,LNLLQ,LNACC,LNREF
300      FORMAT(I2,I1,F12.6,E10.3,E10.3,F5.4,F5.4,F10.4,
     1   F4.2,F8.6,I3,I3,A9,A9,3I1,3I2)
      ELSE
      WRITE(*,*)' PRINTLINE.f :: Invalid data format.Stopping program.'
        STOP
      ENDIF

      RETURN

99    CONTINUE

C=======================================================================
C
C	For debugging purposes ...
C
C=======================================================================
      WRITE(*,*)' '
      WRITE(*,*)' PRINTLINE.f :: Could not print buffer.'
      WRITE(*,*)' BUFFER = '
      WRITE(*,103)BUFFER
103   FORMAT(A)

      print*,LNID,LNISO,LNWAVE,L1,EXP
      print*,LNEINA,EXP1,LNWIDA,LNWIDS
      print*,LNLSE,LNTDEP,LNPSH

      RETURN

      END
