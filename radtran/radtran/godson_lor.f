      real function godson_lor(s,w,v)
C     $Id: godson_lor.f,v 1.2 2011-06-17 15:40:26 irwin Exp $
C     *****************************************************************
C
C     Calculates the equivalent width of a band using Godson Model.
C       Rodgers "Approx. Methods of Calculating Transmission' 1976
C
C     Assumes lines are Lorentz broadened
C       Pat Irwin       18/2/94
C
C     *****************************************************************
      implicit none
      real s,w,v,y,x,pi,i0,i1,a
      real bessi0,bessi1
      parameter (pi=3.1415927)

      x = 2*pi*s*s*v/w
      y = 2*pi*w*w/(s*s)
      i0 = bessi0(y)
      i1 = bessi1(y)

      a = (i0 +2.*y*(i0+i1))*exp(-y) - 1.

      godson_lor = v*x*a

      return
      end
 
      FUNCTION BESSI0(X)
C     $Id: godson_lor.f,v 1.2 2011-06-17 15:40:26 irwin Exp $
      DOUBLE PRECISION Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
     *0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END


      FUNCTION BESSI1(X)
C     $Id: godson_lor.f,v 1.2 2011-06-17 15:40:26 irwin Exp $
      DOUBLE PRECISION Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END

