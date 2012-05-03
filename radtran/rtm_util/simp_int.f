       FUNCTION SIMP_INT(X,NDIM,N,H)
C      **************************************************************
C      Function to integrate array using Simpson's Rule
C
C      Pat Irwin	14/8/00
C
C      **************************************************************
       INTEGER N,I,NDIM
       REAL X(NDIM),SUM,H
       IF(0.5*N.EQ.INT(0.5*N))GOTO 999
       SUM=X(1)+X(N)
       DO 10 I=2,N-1,2
        SUM=SUM+4*X(I)
10     CONTINUE
       DO 20 I=3,N-2,2
        SUM=SUM+2*X(I)
20     CONTINUE
       SIMP_INT=H*SUM/3.0
       RETURN
999    STOP 'Simpson: Integer N must be ODD'
       END
