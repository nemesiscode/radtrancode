      SUBROUTINE noverlapg(idump,delg,ng,ngas,amount,k_gn,dkgndT,
     1 k_g,dkdq)
C***********************************************************************
C_TITL:	NOVERLAPG.f
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
        WRITE(*,*)'NOVERLAPG :: amount = ',(amount(igas),igas=1,ngas)
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
          DO i=1,ng
            k_g1(i) = k_gn(i,igas)
            k_g2(i) = k_gn(i,igas+1)
            dk1dT(i) = dkgndT(i,igas)
            dk2dT(i) = dkgndT(i,igas+1)
          ENDDO

          IF(idump.EQ.1)THEN
            WRITE(*,*)'NOVERLAPG :: igas = ',igas
            WRITE(*,*)'NOVERLAPG :: a1, k1 = ',a1,(k_g1(i),i=1,ng)
            WRITE(*,*)'NOVERLAPG :: a2, k2 = ',a2,(k_g2(i),i=1,ng)
            WRITE(*,*)'NOVERLAPG :: a1, dk1dT = ',a1,(dk1dT(i),i=1,ng)
            WRITE(*,*)'NOVERLAPG :: a2, dk2dT = ',a2,(dk2dT(i),i=1,ng)
          ENDIF

C Skip if first k-distribution = 0.0
          IF(k_g1(ng).EQ.0.0)THEN
            DO I=1,ng
              k_g(I) = k_g2(I)*a2
              dkdq(I,igas) = 0.0
              dkdq(I,igas+1) = k_g2(I)
              dkdq(I,igas+2) = dk2dT(I)*a2
            ENDDO
            GOTO 99
          ENDIF

C Skip if second k-distribution = 0.0
          IF(k_g2(ng).EQ.0.0)THEN
            DO I=1,ng
              k_g(I) = k_g1(I)*a1
              dkdq(I,igas) = k_g1(I)
              dkdq(I,igas+1) = 0.0
              dkdq(I,igas+2) = dk1dT(I)*a1
            ENDDO
            GOTO 99
          ENDIF

          jtest=ng*ng
          if(jtest.gt.maxrank)then
           print*,'Overflow in noverlap.f'
           print*,'Arrays weight and contri will overflow'
           stop
          endif

          nloop = 0
       
          DO I=1,ng
            DO J=1,ng
              nloop = nloop + 1
              weight(nloop) = delg(I)*delg(J)
              contri(nloop) = k_g1(I)*a1 + k_g2(J)*a2
              grad(nloop,igas) = k_g1(I)
              grad(nloop,igas+1) = k_g2(J)
              grad(nloop,igas+2) = dk1dT(I)*a1 + dk2dT(J)*a2
            ENDDO
          ENDDO
        ELSE
C=======================================================================
C
C	Subsequent gases ... add amount*k to previous summed k.
C
C=======================================================================
          a2 = amount(igas+1)
          DO i=1,ng
            k_g1(i) = k_g(i) 
            k_g2(i) = k_gn(i,igas+1)
            dk1dT(i) = dkdq(i,igas+1)           ! dK/dT of previous sum
            dk2dT(i) = dkgndT(i,igas+1)         ! dK/dT of new dist
          ENDDO

          IF(idump.EQ.1)THEN
            WRITE(*,*)'NOVERLAPG :: igas = ',igas
            a1 = 1.0
            WRITE(*,*)'NOVERLAPG :: a1, k1 = ',a1,(k_g1(i),i=1,ng)
            WRITE(*,*)'NOVERLAPG :: a2, k2 = ',a2,(k_g2(i),i=1,ng)
            WRITE(*,*)'NOVERLAPG :: a1, dk1dT = ',a1,(dk1dT(i),i=1,ng)
            WRITE(*,*)'NOVERLAPG :: a2,dk2dT = ',a2,(dk2dT(i),i=1,ng)
          ENDIF

C Skip if first k-distribution = 0.0
          IF(k_g1(ng).EQ.0.0)THEN
            DO I=1,ng
              k_g(I) = k_g2(I)*a2
              DO jgas=1,igas
                dkdq(I,jgas) = 0.0
              ENDDO
              dkdq(I,igas+1) = k_g2(I)
              dkdq(I,igas+2) = dk2dT(I)*a2
            ENDDO
            GOTO 99
          ENDIF

C Skip if second k-distribution = 0.0
          IF(k_g2(ng).EQ.0.0)THEN
            DO I=1,ng
              k_g(I) = k_g1(I)
              dkdq(I,igas+1) = 0.0
              dkdq(I,igas+2) = dk1dT(I)
            ENDDO
            GOTO 99
          ENDIF

          nloop = 0
          DO I=1,ng
            DO J=1,ng
              nloop = nloop + 1
              weight(nloop) = delg(I)*delg(J)
              contri(nloop) = k_g1(I) + k_g2(J)*a2
              DO jgas = 1,igas
                grad(nloop,jgas) = dkdq(I,jgas)
              ENDDO
              grad(nloop,igas+1) = k_g2(J)
              grad(nloop,igas+2) = dk1dT(I) + dk2dT(J)*a2
            ENDDO
          ENDDO
        ENDIF

        CALL rankg(igas+2,delg,ng,nloop,weight,contri,grad,k_g,dkdq)

99      CONTINUE

        IF(idump.EQ.1)THEN
          WRITE(*,*)'NOVERLAPG :: k_g = ',(k_g(i),i=1,ng)
          WRITE(*,*)'NOVERLAPG :: dkdT = ',(dkdq(i,igas+2),i=1,ng)
        ENDIF
100   CONTINUE

      RETURN

      END
