      SUBROUTINE noverlapg1(idump,delg,ng,ngas,amount,k_gn,dkgndT,
     1 k_g,dkdq)
C***********************************************************************
C_TITL:	NOVERLAPG1.f
C
C_DESC:	Combines the absorption coefficient distributions of two
C	overlapping gases. The overlap is implicitly assumed to be random
C	and the k-distributions are assumed to have NG mean values and NG
C	weights. Correspondingly there are NG+1 ordinates in total.
C
C_ARGS:	Input variables:
C	idump			INTEGER	If set to 1, diagnostic info is
C					printed.
C	delg(ng)		REAL	Widths of bins in g-space.
C	ng			INTEGER Number of ordinates.
C	ngas			INTEGER	Number of k-tables to overlap.
C	amount(maxgas)		REAL	Absorber amount of each gas.
C	k_gn(maxg,maxgas)	REAL	K-distributions of the different
C					gases.
C	dkgndT(maxg,maxgas)	REAL	dK/dTemperature for the different
C					gases.
C
C	Output variables:
C	k_g(maxg)		REAL	K-distribution of combined gases.
C	dkdq(maxg,maxgas+1)	REAL	dk/dparam: 1 --> ngas is dk/dvmr,
C					ngas+1 is dk/dTemp.
C
C_FILE:	No files openned.
C
C_CALL:	rankg	Sorts randomised k-coefficients, and associated
C		gradients into the mean k-distribution and gradient.
C
C_HIST:	15/2/95	PGJI	Original version.
C	27/7/01	PGJI	Updated for gradients.
C	29/2/12	PGJI	Updated for Radtrans2.0
C
C***************************** VARIABLES *******************************
 
      IMPLICIT NONE

C The include file ...
      INCLUDE '../includes/arrdef.f'
C ..//includes/arrdefs.f defines the maximum values for a series of
C variables (layers, bins, paths, etc.)


C The input variables ...
      INTEGER idump,ng,ngas,jtest
      REAL delg(maxg),k_g(maxg),amount(maxgas)
      REAL k_gn(maxg,maxgas),dkgndT(maxg,maxgas),dkdq(maxg,maxgas+1)


C General variables ...
      INTEGER i,j,nloop,igas,jgas
      REAL dk1dT(maxg),dk2dT(maxg),k_g1(maxg),k_g2(maxg)
      REAL contri(maxrank),weight(maxrank),grad(maxrank,maxgas+1),a1,a2

C******************************** CODE *********************************

      IF(idump.EQ.1)THEN
        WRITE(*,*)'NOVERLAPG1 :: amount = ',(amount(igas),igas=1,ngas)
      ENDIF

      IF(NG.NE.1)THEN
       print*,'NOVERLAPG1 - NG Must be equal to 1'
       stop
      ENDIF

      DO 100 igas=1,ngas-1
        IF(igas.EQ.1)THEN
C=======================================================================
C
C	First pair of gases
C
C=======================================================================
          a1 = amount(igas)
          a2 = amount(igas+1)
          k_g1(1) = k_gn(1,igas)
          k_g2(1) = k_gn(1,igas+1)
          dk1dT(1) = dkgndT(1,igas)
          dk2dT(1) = dkgndT(1,igas+1)


C Skip if first k-distribution = 0.0
          IF(k_g1(1).EQ.0.0)THEN
              k_g(1) = k_g2(1)*a2
              dkdq(1,igas) = 0.0
              dkdq(1,igas+1) = k_g2(I)
              dkdq(1,igas+2) = dk2dT(I)*a2
            GOTO 99
          ENDIF

C Skip if second k-distribution = 0.0
          IF(k_g2(1).EQ.0.0)THEN
              k_g(1) = k_g1(1)*a1
              dkdq(1,igas) = k_g1(1)
              dkdq(1,igas+1) = 0.0
              dkdq(1,igas+2) = dk1dT(1)*a1
            GOTO 99
          ENDIF


          nloop = 1
       
          k_g(1) = k_g1(1)*a1 + k_g2(1)*a2
          dkdq(1,igas) = k_g1(1)
          dkdq(1,igas+1) = k_g2(1)
          grad(1,igas+2) = dk1dT(1)*a1 + dk2dT(1)*a2

        ELSE
C=======================================================================
C
C	Subsequent gases ... add amount*k to previous summed k.
C
C=======================================================================
          a2 = amount(igas+1)
          k_g1(1) = k_g(1) 
          k_g2(1) = k_gn(1,igas+1)
          dk1dT(2) = dkdq(1,igas+1)           ! dK/dT of previous sum
          dk2dT(2) = dkgndT(1,igas+1)         ! dK/dT of new dist

C Skip if first k-distribution = 0.0
          IF(k_g1(1).EQ.0.0)THEN
              k_g(1) = k_g2(1)*a2
              DO jgas=1,igas
                dkdq(1,jgas) = 0.0
              ENDDO
              dkdq(1,igas+1) = k_g2(1)
              dkdq(1,igas+2) = dk2dT(1)*a2
            GOTO 99
          ENDIF

C Skip if second k-distribution = 0.0
          IF(k_g2(ng).EQ.0.0)THEN
              k_g(1) = k_g1(1)
              dkdq(1,igas+1) = 0.0
              dkdq(1,igas+2) = dk1dT(1)
            GOTO 99
          ENDIF

          k_g(1) = k_g1(1) + k_g2(1)*a2
          DO jgas = 1,igas
                dkdq(1,jgas) = dkdq(1,jgas)
          ENDDO
          dkdq(1,igas+1) = k_g2(1)
          dkdq(1,igas+2) = dk1dT(1) + dk2dT(1)*a2

        ENDIF

99      CONTINUE

        IF(idump.EQ.1)THEN
          WRITE(*,*)'NOVERLAPG1 :: k_g = ',(k_g(i),i=1,ng)
          WRITE(*,*)'NOVERLAPG1 :: dkdT = ',(dkdq(i,igas+2),i=1,ng)
        ENDIF
100   CONTINUE

      RETURN

      END
