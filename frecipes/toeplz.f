      SUBROUTINE TOEPLZ(R,X,Y,N)
      PARAMETER (NMAX=100)
      DIMENSION R(*),X(N),Y(N),G(NMAX),H(NMAX)
      IF(R(N).EQ.0.) GO TO 99
      X(1)=Y(1)/R(N)
      IF(N.EQ.1)RETURN
      G(1)=R(N-1)/R(N)
      H(1)=R(N+1)/R(N)
      DO 15 M=1,N
        M1=M+1
        SXN=-Y(M1)
        SD=-R(N)
        DO 11 J=1,M
          SXN=SXN+R(N+M1-J)*X(J)
          SD=SD+R(N+M1-J)*G(M-J+1)
11      CONTINUE
        IF(SD.EQ.0.)GO TO 99
        X(M1)=SXN/SD
        DO 12 J=1,M
          X(J)=X(J)-X(M1)*G(M-J+1)
12      CONTINUE
        IF(M1.EQ.N)RETURN
        SGN=-R(N-M1)
        SHN=-R(N+M1)
        SGD=-R(N)
        DO 13 J=1,M
          SGN=SGN+R(N+J-M1)*G(J)
          SHN=SHN+R(N+M1-J)*H(J)
          SGD=SGD+R(N+J-M1)*H(M-J+1)
13      CONTINUE
        IF(SD.EQ.0..OR.SGD.EQ.0.)GO TO 99
        G(M1)=SGN/SGD
        H(M1)=SHN/SD
        K=M
        M2=(M+1)/2
        PP=G(M1)
        QQ=H(M1)
        DO 14 J=1,M2
          PT1=G(J)
          PT2=G(K)
          QT1=H(J)
          QT2=H(K)
          G(J)=PT1-PP*QT2
          G(K)=PT2-PP*QT1
          H(J)=QT1-QQ*PT2
          H(K)=QT2-QQ*PT1
          K=K-1
14      CONTINUE
15    CONTINUE
      PAUSE 'never get here'
99    PAUSE 'Levinson method fails: singular principal minor'
      END
