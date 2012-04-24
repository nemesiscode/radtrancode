      SUBROUTINE REMSP(TEXT)
C     $Id: remsp.f,v 1.1 2000-08-17 09:26:54 irwin Exp $
C     **********************************************************************
C     Removes leading spaces from text string
C
C	1/1/90	SBC	Original Version
C	3/10/94	PGJI	Added $Id: remsp.f,v 1.1 2000-08-17 09:26:54 irwin Exp $
C
C     **********************************************************************
      CHARACTER TEXT*(*)
      INTEGER I,J,K,L
      J=LEN(TEXT)
      DO 100 I=1,J
      IF(TEXT(I:I).NE.' ')GOTO 10
100   CONTINUE
      RETURN
10    CONTINUE
      IF(I.EQ.1)RETURN
      DO 20 K=I,J
      L=K-I+1
      TEXT(L:L)=TEXT(K:K)
20    CONTINUE
      DO 30 K=J-I+2,J
      TEXT(K:K)=' '
30    CONTINUE
      END
