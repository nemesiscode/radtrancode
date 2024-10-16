C     $Id: dbcom.f,v 1.6 2011-06-17 14:49:05 irwin Exp $
C***********************************************************************
C_TITL: DBCOM
C
C_DESC:	Common variables used by all linedata routines. Mixed data types
C	are avoided.
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************

      INTEGER MAXDGAS,MAXISO
      PARAMETER (MAXDGAS=150,MAXISO=20)

      INTEGER*8 DBSIZ
      INTEGER DBRECL,DBFORM,LOCID(MAXDGAS)
      INTEGER DBISO(MAXISO,MAXDGAS),DBNISO(MAXDGAS),DBREC,DBLUN,INDRL
      COMMON /DBINT/ DBSIZ,DBRECL,DBFORM,LOCID,DBISO,DBNISO,
     1 DBREC,DBLUN,INDRL

      REAL MASSNO(MAXISO,MAXDGAS),RELABU(MAXISO,MAXDGAS)
      REAL DBQTA(MAXISO,MAXDGAS),DBQTB(MAXISO,MAXDGAS)
      REAL DBQTC(MAXISO,MAXDGAS),DBQTD(MAXISO,MAXDGAS)
      COMMON /DBREAL/ MASSNO,RELABU,DBQTA,DBQTB,DBQTC,DBQTD

      CHARACTER*16 DBNAME
      CHARACTER*8 GASNAM(MAXDGAS)
      CHARACTER*1 LNFLAG
      CHARACTER*128 IPFILE,DBFILE,KEYFIL,INDFIL,ISOFIL,GASFIL,XSCFIL
      CHARACTER*128 DGFILE,RADFILE
      COMMON /DBCHAR/ IPFILE,DBFILE,KEYFIL,INDFIL,ISOFIL,GASFIL,XSCFIL,
     1 GASNAM,DBNAME,DGFILE,RADFILE

      INTEGER LNID,LNISO,LNUGQI,LNLGQI,LNACC(3),LNREF(3)
      INTEGER LNACC04(6),LNREF04(6),LNUJ,LNLJ,LNUK,LNLK
      COMMON /LNINT/ LNID,LNISO,LNUGQI,LNLGQI,LNACC,LNREF,
     1  LNACC04,LNREF04,LNUJ,LNLJ,LNUK,LNLK

      DOUBLE PRECISION LNWAVE,LNSTR,LNPROB,LNEINA
      REAL LNWIDA,LNWIDS,LNLSE,LNTDEP,LNTDEPS
      REAL LNPSH,LDOUBV,UWGHT,LWGHT,LNTDEP1
      REAL LNWIDA1

      COMMON /LNREAL/LNWIDA,LNWIDA1,LNWIDS,LNLSE,LNTDEP,LNTDEP1,
     1LNTDEPS,LNPSH,LDOUBV,UWGHT,LWGHT
      COMMON /LNDBLE/LNWAVE,LNSTR,LNPROB,LNEINA

      CHARACTER LNULQ*9,LNLLQ*9,TRS1*9,TRS2*9,RN1*9,RN2*9
      CHARACTER LNUGQI04*15,LNLGQI04*15,LNULQ04*15,LNLLQ04*15
      CHARACTER*18 LNUVQ,LNLVQ
      COMMON /LNCHAR/LNULQ,LNLLQ,TRS1,TRS2,RN1,RN2,
     1 LNUGQI04,LNLGQI04,LNULQ04,LNLLQ04,LNFLAG,LNUVQ,LNLVQ
C--------------------------------------------------------------
