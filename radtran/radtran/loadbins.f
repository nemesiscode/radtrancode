      SUBROUTINE LOADBINS(WING,NGAS,IDGAS,ISOGAS)
C     $Id: loadbins.f,v 1.5 2011-06-17 15:40:27 irwin Exp $
C***********************************************************************
C_TITL:	LOADBINS.f
C
C_DESC:	Loads linedata "bins" for GENLBL.
C
C_ARGS:	Input variables:
C	WING	REAL	Bin size.
C	NGAS	INT	Number of gases.
C	IDGAS	INT	Local gas identifier array.
C	ISOGAS	INT	Local isotope identifier array.
C
C_FILE:	No files openned.
C
C_CALL:	liness	Reads in linedata bins for program LBL.
C
C_HIST:	7oct92	SBC	ORIGINAL VERSION added this header to existing
C			routine.
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
C ../includes/dbcom.f stores the line database variables.
      INCLUDE '../includes/dbcom.f'
C ../includes/bincom.f stores the line bin variables (including NLINES,
C FSTLIN, LSTLIN) and band parameters.
      INCLUDE '../includes/bincom.f'
C ../includes/lincom.f stores the linedata variables (MAXLIN, VLIN, SLIN,
C ALIN, ELIN, SBLIN, TDW, TDWS and that lot).
      INCLUDE '../includes/lincom.f'

      INTEGER I,J
      INTEGER NLIN,FIRST,LAST
      INTEGER NGAS,IDGAS(NGAS),ISOGAS(NGAS)
      DOUBLE PRECISION VMIN,VMAX,DELV
      REAL WING

      REAL TINW(MAXDGAS)
      DATA (TINW(J),J=1,42)/0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5/

      LOGICAL WARNTD

C******************************** CODE *********************************

      WARNTD= .TRUE.

      VMIN= DBLE(VBIN(1))
      VMAX= DBLE(VBIN(NBIN) + WING)
      DELV= VMAX - VMIN

      FIRST= 1
      LAST= 0
      NLIN= 0

      WRITE(*,*)' LOADBINS.f :: nbin, wing = ',nbin,wing
      WRITE(*,*)' LOADBINS.f :: vmin, vmax, delv = ',vmin,vmax,delv
      WRITE(*,*)' LOADBINS.f :: maxlin = ',maxlin
      WRITE(*,*)' LOADBINS.f :: ngas = ',ngas
      DO i=1,ngas
        WRITE(*,*)' LOADBINS.f :: idgas, isogas = ',idgas,isogas
      ENDDO

      CALL LINESS(VMIN,DELV,MAXLIN,VLIN,SLIN,ALIN,ELIN,IDLIN,
     1 SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,FIRST,LAST,NGAS,
     2 IDGAS,ISOGAS)

      WRITE(*,*)' LOADBINS.f :: Number of lines = ',NLIN
      WRITE(*,*)' LOADBINS.f :: First, last line= ',FIRST,LAST

      DO 11 I=1,NBIN
        NLINES(I)= 0
        FSTLIN(I)= 0
        LSTLIN(I)= 0
11    CONTINUE

      DO 51 J=1,NLIN
        IF(TDW(J).LT.1.E-30)THEN
          TDW(J)= TINW(IDLIN(J))
          IF(WARNTD)THEN
            WRITE(*,*)' LOADBINS.f :: Warning: setting default line'
            WRITE(*,*)' width temperature exponents.'
            WARNTD= .FALSE.
          ENDIF
        ENDIF

        IF(TDWS(J).LT.1.E-30)THEN
          TDWS(J)= TINW(IDLIN(J))
          IF(WARNTD)THEN
            WRITE(*,*)' LOADBINS.f :: Warning: setting default line'
            WRITE(*,*)' width temperature exponents.'
            WARNTD= .FALSE.
          ENDIF
        ENDIF

        I= 1 + INT((VLIN(J) - VBIN(1))/WING)
        NLINES(I)= NLINES(I) + 1
        IF(FSTLIN(I).EQ.0)FSTLIN(I)= J
        LSTLIN(I)= J
51    CONTINUE

      DO 25 I=1,NBIN
        IF(FSTLIN(I).EQ.0)FSTLIN(I)= 1
25    CONTINUE

      PRINT*,'LOADBINS: NBIN = ',NBIN
      PRINT*,'I,VV,NLINES(I),FSTLIN(I),LSTLIN(I)'
      DO 26 I=1,NBIN
        WRITE(*,*)I,VMIN+(I-1)*WING,NLINES(I),FSTLIN(I),LSTLIN(I)
26    CONTINUE

      RETURN

      END
