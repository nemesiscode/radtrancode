C     **********************************************************
C     Parameters defining new continuum bins and also layer temperature
C     properties.

      INTEGER NBIN
      REAL VBIN(MAXBIN)
      REAL CONTINK(IORDP1,MAXLAY,MAXBIN)   
      REAL CONVALS(IORDP1,MAXLAY,MAXBIN),CONWAV(IORDP1)
      REAL MATRIX(IORDP1,IORDP1),UNIT(IORDP1,IORDP1)
      REAL TCORS1(MAXLAY,MAXGAS),TCORS2(MAXLAY)
      REAL TCORDW(MAXLAY,MAXGAS)
C
      INTEGER FSTLIN(2,MAXBIN),LSTLIN(2,MAXBIN),LASTBIN(2)

C     IORDER is the order of the continuum polynomial.
C     IORDP1 is the number of parameters in polynomial fit.
C     CONTINK holds continuum polynomial for each bin.
C     CONVALS holds the absorption coefficients at wavenumbers CONWAV prior to
C     fitting continuum.
C     MATRIX and UNIT are both used for matrix inversion when computing
C     polynomials where insufficient points for least squares.

      COMMON /CONCOMB/NBIN,VBIN,CONTINK,CONVALS,CONWAV,MATRIX,UNIT,
     1 FSTLIN,LSTLIN,LASTBIN,TCORS1,TCORS2,TCORDW

