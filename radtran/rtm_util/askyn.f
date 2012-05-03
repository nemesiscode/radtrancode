      LOGICAL FUNCTION ASKYN(TEXT)
C     $Id: askyn.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ********************************************************************
C     Passes a prompt to user and returns true for a Y or y input
C     or false otherwise
C
C     1/1/90	Original Version:	SBC
C     3/10/94	Updated Header		PGJI
C
C     ********************************************************************
      CHARACTER TEXT*(*)
      CHARACTER YN

      CALL PROMPT(TEXT)
      READ(*,10)YN
10    FORMAT(1A1)
      CALL UPCASE(YN)
      IF(YN.EQ.'Y')THEN
        ASKYN=.TRUE.
       ELSE
        ASKYN=.FALSE.
        END IF
      RETURN
      END
