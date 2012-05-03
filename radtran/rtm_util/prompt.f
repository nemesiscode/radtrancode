      SUBROUTINE PROMPT(CHAR)
C     $Id: prompt.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C---------------------------------------------------------------------------
C
C_TITLE: PROMPT: simply sends prompts to terminal
C
C_ARGS:  None.
C
C_KEYS:  PROG, VMS, CAL .
C
C_DESCR: sends prompt to terminal. Used to aid transportability
C        $ format used is not standard FORTRAN 77 so should not
C        be included in transportable code
C
C_FILES: 
C
C_CALLS: 
C                          
C_BUGS:  
C
C_HIST:   ORIGINAL VERSION Mark Radford for NIMS Calibration Software
C         29jan90 SBC mod. to ignore trailing spaces
C	  03oct94 PGJI Added $Id: prompt.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C
C_END:
C
C--------------------------------------------------------------------------
      INTEGER K
      CHARACTER*(*) CHAR
      DO 100 K=LEN(CHAR),1,-1
      IF(CHAR(K:K).NE.' ')THEN
        IF(CHAR(K:K).EQ.':')THEN
          WRITE(*,10) CHAR(1:K)
10      FORMAT(' ',A,$)
         ELSE
C         forcing prompt character to :
          WRITE(*,11) CHAR(1:K)
11        FORMAT(' ',A,' :',$)
          END IF
        RETURN
       END IF
100   CONTINUE
      RETURN
      END
