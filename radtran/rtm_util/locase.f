      SUBROUTINE LOCASE(Q)
C     $Id: locase.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ********************************************************************
C     Forces the character string Q to lower case
C
C     1/1/90    Original Version:       SBC
C     3/10/94   Updated Header          PGJI
C
C     ********************************************************************
      CHARACTER Q*(*)
      INTEGER I,J
      DO 100 I=1,LEN(Q)
      J=ICHAR(Q(I:I))
      IF(J.GE.65.AND.J.LE.90)J=J+32
      Q(I:I)=CHAR(J)
100   CONTINUE
      RETURN
      END
