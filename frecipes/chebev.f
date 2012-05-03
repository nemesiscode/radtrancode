      FUNCTION CHEBEV(A,B,C,M,X)
      DIMENSION C(M)
      IF ((X-A)*(X-B).GT.0.) PAUSE 'X not in range.'
      D=0.
      DD=0.
      Y=(2.*X-A-B)/(B-A)
      Y2=2.*Y
      DO 11 J=M,2,-1
        SV=D
        D=Y2*D-DD+C(J)
        DD=SV
11    CONTINUE
      CHEBEV=Y*D-DD+0.5*C(1)
      RETURN
      END
