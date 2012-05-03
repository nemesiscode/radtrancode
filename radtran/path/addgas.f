      SUBROUTINE ADDGAS(ID,ISO,LOCID)
C     $Id: addgas.f,v 1.1.1.1 2000-08-17 09:26:56 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  title: 
C
C_KEYS:   
C
C_DESCR:  
C
C_ARGS:   
C
C_FILES : unit n - description
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   23jun92 SBC Original version copied out of PATH
C
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C     local variables
      INTEGER J,ID,ISO,LOCID,ILAY
C     search gas arrays to see if gas already included
      DO 405 J=1,NGAS
      IF(IDGAS(J).EQ.ID.AND.ISOGAS(J).EQ.ISO)THEN
        LOCID=J
        GOTO 406
        END IF
405   CONTINUE
C     if can"t fing gas then insert it into gas arrays
      NGAS=NGAS+1
      LOCID=NGAS
      IDGAS(NGAS)=ID
      ISOGAS(NGAS)=ISO
C     if it"s a new gas then its amount must be set to zero in the old
C     layers
      DO 407 ILAY=1,NLAYER-1
      AMOUNT(ILAY,NGAS)=0.
      PP(ILAY,NGAS)=0.
407   CONTINUE
406   CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------------

