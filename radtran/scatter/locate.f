      SUBROUTINE LOCATE(XX,N,X,J)
      DIMENSION XX(N)
C     Given an array xx(1:n), and given a value x, returns a value
C     j such that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, 
C     either increasing or decreasing. j=0 or j=n is returned to indicate that
C     x is out of range.
    
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GE.XX(1)).EQV.(X.GE.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      IF(X.EQ.XX(1))THEN
       J=1
      ELSE IF(X.EQ.XX(N))THEN
       J=N-1
      ELSE
       J=JL
      ENDIF

      RETURN
      END
