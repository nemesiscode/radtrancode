      SUBROUTINE QUAD3D(X1,X2,SS)
      EXTERNAL H
      CALL QGAUSX(H,X1,X2,SS)
      RETURN
      END
C     
      FUNCTION F(ZZ)
      EXTERNAL FUNC
      COMMON /XYZ/ X,Y,Z
      Z=ZZ
      F=FUNC(X,Y,Z)
      RETURN
      END
C     
      FUNCTION G(YY)
      EXTERNAL Z1,Z2,F
      COMMON /XYZ/ X,Y,Z
      Y=YY
      CALL QGAUSZ(F,Z1(X,Y),Z2(X,Y),SS)
      G=SS
      RETURN
      END
C     
      FUNCTION H(XX)
      EXTERNAL Y1,Y2,G
      COMMON /XYZ/ X,Y,Z
      X=XX
      CALL QGAUSY(G,Y1(X),Y2(X),SS)
      H=SS
      RETURN
      END
