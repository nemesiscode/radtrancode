	SUBROUTINE LOBATTO( X, W, NMU)
C       COMPUTE LOBATTO QUADRATURE ABSCISSAS AND WEIGHTS.
C       ABSCISSAS ARE POSITIVE ZEROES OF LEGENDRE POLYNOMIAL DERIVATIVES.
	PARAMETER (MAX=45)
	DOUBLE PRECISION X(nmu), W(nmu)
	DOUBLE PRECISION A(MAX), REZ(MAX), IMZ(MAX), TOL, SUM, POLY
	DOUBLE PRECISION icoef((max+1)/2),jcoef(max/2+1)
	REAL*4 DUM/1./
	integer count

        nn=nmu*2
	N = NN-1

C       X(I) IS THE (I-1)ST ZERO OF P'(NN-1).
	IF (N.LT.2 .OR. N.GT.MAX-1) GO TO 900
	call legendre(n,icoef,jcoef)

C       INITIALIZE:
	DO I=1,N
	  A(I) = 0.0
	ENDDO

	NT = (N+1)/2		! # TERMS IN THE POLYNOMIAL
	N0 = 2*(NT-1)		! DEGREE OF THE POLYNOMIAL
	N1 = N0+1		! FOR C02AEF
	DO I=1,NT
	  A(2*I-1) = ICOEF(NT+1-I)
	ENDDO
	TOL = X02AAF( DUM)
	IFAIL = 0
	CALL C02AEF( A, N1, REZ, IMZ, TOL, IFAIL)
	IF (IFAIL.EQ.1) GO TO 901
	ITER = 1
	DO WHILE (IFAIL.EQ.2)
	  ITER = ITER+1
	  TOL = X02AAF( DUM)
	  CALL C02AEF( A, N1, REZ, IMZ, TOL, IFAIL)
	  IF (IFAIL.EQ.1) GO TO 901
	  IF (ITER.GT.45) GO TO 902
	ENDDO
	X((n+1)/2) = 1.0
	W(n/2+1) = 2./(DFLOAT(N*(N+1)))
	NN = 1
	DO I=2,N0+1
	  IF (REZ(I-1).GT.0.) THEN	! DISCARD THE SYMMETRIC NEG. ROOTS
	    NN = NN+1
	    nx=(n+3)/2-nn
	    if (nx.gt.0) X(nx) = REZ(I-1)
	    POLY = 1.0
	    IF (N/2.EQ.(N-1)/2.and.nx.gt.0) POLY = X(nx)
	    SUM = JCOEF(1)*POLY
	    DO J = 1,N/2
	      if (nx.gt.0) then
	        POLY = POLY*X(nx)*X(nx)
	      else
	        poly = 0.
	      endif
	      SUM = SUM+JCOEF(J+1)*POLY
	    ENDDO
	    W(n/2+2-nn) = 2./(DFLOAT(N*(N+1))*SUM*SUM)
	  ENDIF
	ENDDO
	IF (N0+1.EQ.N-1) THEN
	  NN = NN+1
	  X(NN) = 0.0
	  SUM = 1./JCOEF(1)
	  W(n/2+2-nn) = SUM*SUM/DFLOAT(N*(N+1))
	ENDIF

c  Sort into order
	count=1
	do while (count.ge.1)
	  count=0
	  do i=1,nmu-1
	    if (x(i).gt.x(i+1)) then
	      dum = x(i)
	      x(i) = x(i+1)
	      x(i+1) = dum
	      dum = w(i)
	      w(i) = w(i+1)
	      w(i+1) = dum
	      count=1
	    endif
	  enddo
	enddo 

	RETURN

900	PRINT*,' N IS TOO LARGE, MAX =', MAX
	GO TO 999
901	PRINT*,' ERROR, IFAIL=1'
	GO TO 999
902	PRINT*,' MORE THAN 45 ITERATIONS REQUIRED'

999	CALL EXIT
	END

	subroutine legendre(nn,icoeff,jcoeff)
c
c  Calculates legendre polynomial coefficients.
c
	implicit none
	integer k,n,nn,nmax,m,mmax,max
	parameter (max=45)
	real*8 icoeff((max+1)/2),jcoeff(max/2+1)
	double precision rcoeff(0:(max+1)/2)
	double precision fact

          n=nn
	  mmax=n/2
	  do m=0,mmax
	    rcoeff(m)=2.**n*fact(m)*fact(n-m)*fact(n-2*m)
	    rcoeff(m)=(-1.)**m*fact(2*n-2*m)/rcoeff(m)
	  enddo
	  do m=0,mmax
	    jcoeff(mmax+1-m)=dble(rcoeff(m))
	  enddo
	  if(mmax*2.ne.n)then
	    do m=0,mmax+1
	      icoeff(m+1)=jcoeff(m+1)*dble(1+(2*m))
	    enddo
	  else
	    do m=0,mmax
	      icoeff(m)=jcoeff(m+1)*dble(2*m)
	    enddo
	  endif
	return
  
        end
              
	double precision function fact(m)      
c
c Factorial worked out. Strange functions are for conversions
c to and from REAL*16.
c
	implicit none
	integer k,m

	fact=1.
	if (m.eq.0) goto 1
	do k=1,m
	  fact=fact*dble(floatj(k))
	enddo
1	end	
