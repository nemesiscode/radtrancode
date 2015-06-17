      SUBROUTINE rankg(nparam,gw,ng,nloop,weight,cont,grad,k_g,dkdq)
C***********************************************************************
C_TILE:	RANKG.f
C
C_DESC:	Subroutine to sort randomised k-coefficients, and associated
C	gradients into the mean k-distribution and gradient.
C	
C_ARGS:	Input variables:
C	nparam			INTEGER	Number of gradients to consider
C	gw(maxg)		REAL	Required weights of final k-dist.
C	ng			INTEGER	Number of weights.
C	nloop			INTEGER	Number of points in randomised
C					k-distribution.
C	weight(maxrank)		REAL	Weights of points in random k-dist
C	cont(maxrank)		REAL	Random k-coeffs.
C	grad(maxrank,maxgas+1)	REAL	Gradients of random k-coeffs.
C
C	Output variables
C	k_g(maxg)		REAL	Mean k-dist.
C	dkdq(maxg)		REAL	Mean gradients.
C
C_FILE:	No files openned.
C
C_CALL:	sort2g	Sorts the input vector into ascending order.
C
C_HIST:	23/2/96	Pat Irwin	Original
C	2/3/96	A.L.Weir
C	9/5/01	Pat Irwin
C	31/7/01	Pat Irwin	Commented
C	29/2/12	Pat Irwin	Updated for Radtrans2.0
C***************************** VARIABLES *******************************

      IMPLICIT NONE

C The include file ...
      INCLUDE '../includes/arrdef.f'
C ../includes/arrdef.f defines the maximum values for a series of
C variables (layers, bins, paths, etc.)


C The input and output variables ...
      INTEGER nparam,ng,nloop
      REAL gw(maxg),k_g(maxg),dkdq(maxg,maxgas+1)
      REAL grad(maxrank,maxgas+1),cont(maxrank),weight(maxrank)
      REAL tmp(maxrank)

C General variables ...
      INTEGER i,ig,iparam,ico(maxrank)
      DOUBLE PRECISION g_ord(maxg), gdist(0:maxrank), sum, frac

C******************************** CODE *********************************

C=======================================================================
C
C	Sum delta gs to get cumulative g ordinate. rank cont and weight
C	arrays in ascending order of k (i.e. cont values).
C
C=======================================================================
      g_ord(1) = 0.0
      DO I=1,ng
        g_ord(I+1) = g_ord(I) + gw(I)
      ENDDO
C     Make sure g_ord(ng+1)=1. (rounding errors can lead to numbers just
C                               less than 1.0)
      if(g_ord(ng+1).lt.1.0)g_ord(ng+1)=1.


      DO I=1,nloop
        ico(I) = I
      ENDDO

C Sort random k-coeffs into order. Integer array ico records which swaps
C have been made so that we can also re-order the gradients.
      CALL sort2g(nloop, cont, ico)

C Resort the weights:
      DO I=1,nloop
        tmp(I) = weight(ico(I))
      ENDDO
      DO I=1,nloop
        weight(I) = tmp(I)
      ENDDO

C Resort the gradients
      DO iparam=1,nparam
        DO I=1,nloop
          tmp(I) = grad(ico(I),iparam)
        ENDDO
        DO I=1,nloop
          grad(I,iparam) = tmp(I)
        ENDDO
      ENDDO 


C=======================================================================
C
C	Now form new g(k) and gradients by summing over weight. The new
C       k(g) is then found by ascending the arranged ks and getting the
C       weighted averages of the k values in the relevant g
C       interval. Linear interpolation is used to deal with steps crossing
C       over the required g intervals.
C
C=======================================================================
      gdist(0) = 0.0
      gdist(1) = weight(1)
      DO I=2,nloop
        gdist(I) = weight(I) + gdist(I-1) 
      ENDDO

      DO I=1,ng
        k_g(I) = 0.0
        DO iparam=1,nparam
          dkdq(I,iparam) = 0.0
        ENDDO
      ENDDO

      ig = 1
      sum = 0.0
      DO I=1,nloop
        IF(gdist(I).LT.g_ord(ig+1))THEN
          k_g(ig) = k_g(ig) + cont(I) * weight(I)
          DO iparam=1,nparam
            dkdq(ig,iparam) = dkdq(ig,iparam) + grad(I,iparam)*weight(I)
          ENDDO
          sum = sum + weight(I)
        ELSE
          frac = (g_ord(ig+1) - gdist(I-1))/(gdist(I) - gdist(I-1))
          k_g(ig) = k_g(ig) + sngl(frac)*cont(I)*weight(I)
          DO iparam=1,nparam
            dkdq(ig,iparam) = dkdq(ig,iparam) + 
     1      sngl(frac)*grad(I,iparam)*weight(I)
          ENDDO
          sum = sum + frac*weight(I)
          k_g(ig) = k_g(ig)/sngl(sum)
          DO iparam=1,nparam
            dkdq(ig,iparam) = dkdq(ig,iparam)/sngl(sum)
          ENDDO
          ig = ig + 1
          sum = (1.0 - frac)*weight(I)
          k_g(ig) = sngl(1.0 - frac)*cont(I)*weight(I)
          DO iparam=1,nparam
            dkdq(ig,iparam) = sngl(1.-frac)*grad(I,iparam)*weight(I)
          ENDDO
        ENDIF
      ENDDO

      IF(ig.EQ.ng)THEN
         k_g(ig) = k_g(ig)/sngl(sum)
         DO iparam=1,nparam
           dkdq(ig,iparam) = dkdq(ig,iparam)/sngl(sum)
         ENDDO
      ENDIF

      RETURN

      END
