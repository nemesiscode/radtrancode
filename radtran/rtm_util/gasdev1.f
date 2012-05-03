      REAL FUNCTION GASDEV1(IDUM)
C     ***************************************************************
C     Returns a randum number with a normal gaussian distribution of unit
C     standard deviation.
C 
C     Pat Irwin		14/8/00
C 
C     ***************************************************************
      INTEGER ISET
      REAL GSET

      COMMON /GSTORE/ISET,GSET

      IF(ISET.NE.0)THEN     
1       V1=2.*RAN11(IDUM)-1.
        V2=2.*RAN11(IDUM)-1.
        R=V1**2 + V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV1=V2*FAC
        ISET=0
      ELSE
        GASDEV1=GSET
        ISET=1
      ENDIF

      RETURN
      END
