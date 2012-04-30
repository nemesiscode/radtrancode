      SUBROUTINE PEARSN(X,Y,N,R,PROB,Z)
      PARAMETER (TINY=1.E-20)
      DIMENSION X(N),Y(N)
      AX=0.
      AY=0.
      DO 11 J=1,N
        AX=AX+X(J)
        AY=AY+Y(J)
11    CONTINUE
      AX=AX/N
      AY=AY/N
      SXX=0.
      SYY=0.
      SXY=0.
      DO 12 J=1,N
        XT=X(J)-AX
        YT=Y(J)-AY
        SXX=SXX+XT**2
        SYY=SYY+YT**2
        SXY=SXY+XT*YT
12    CONTINUE
      R=SXY/SQRT(SXX*SYY)
C      Z=0.5*ALOG(((1.+R)+TINY)/((1.-R)+TINY))
      Z=1.0
      DF=N-2
      T=R*SQRT(DF/(((1.-R)+TINY)*((1.+R)+TINY)))
      PROB=BETAI(0.5*DF,0.5,DF/(DF+T**2))
C     PROB=ERFCC(ABS(Z*SQRT(N-1.))/1.4142136)
      RETURN
      END
