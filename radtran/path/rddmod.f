      SUBROUTINE RDDMOD(TEXT)
C     $Id: rddmod.f,v 1.5 2004-11-25 11:53:17 irwin Exp $
C***********************************************************************
C_TITL:	RDDMOD
C
C_DESC:	Reads in a dust profiles for path.f.
C
C_ARGS:	Input variable:
C	text	CHARACTER*(*)	Text string containing the name of the
C				profile runname.prf file.
C				CHARACTER*(*) declares an incoming
C				character variable whose length is
C				unknown.
C
C_FILE:	unit=1	Atmospheric model (.prf)
C
C_CALL:	remsp		Remove spaces from input string.
C       file		Forces file extension.
C       locase		Make lower case the input string.
C	verint		Performs vertical interpolation and integration
C			of atmospheric profiles
C
C_HIST:	3Dec93	PGJI	Original version.
C***********************************************************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      include '../includes/emcee.f'
C ../includes/pathcom.f holds the variables used by the software when
C calculating atmospheric paths (e.g. NPATH and IMOD).
      INCLUDE '../includes/laycom.f'
C ../includes/laycom.f holds variables used only by the path software
C parameters are passed between routines mostly using common blocks
C because of the extensive use of large arrays. NOTE: laycom uses
C parameters defined in pathcom.

      INTEGER i,j,k,l,nn,ilun
C N: Maximum VMR which can be read in from the next block.
C ILUN: File unit used for openning files.

      REAL TD1,TDUST(MAXPRO)

      CHARACTER*(*) text
      CHARACTER*100 buffer,ipfile
C IPFILE: Input filename.

C********************************* CODE ********************************

C Reading in in vertical profiles produced by dust_profile.f
      ilun = 1
      READ(text,1)ipfile
1     FORMAT(A)
      CALL remsp(ipfile)
      CALL locase(ipfile) 
      CALL file(ipfile,ipfile,'prf')
      WRITE(*,*)' RDDMOD.f :: reading dust-model: ',ipfile
      OPEN(UNIT=ilun,FILE=ipfile,STATUS='OLD')
C First skip header (if any)
54    READ(ILUN,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)NN,NDUST
      WRITE(*,*)' RDDMOD.f :: dust-model has ',ndust,' aerosol types'
      NCONT = NDUST
      
      DO 105 J=1,NN
         READ(ILUN,*)DUSTH(J),(DUST(I,J),I=1,NDUST)
C        WRITE(*,*)DUSTH(J),(DUST(I,J),I=1,NDUST)
105   CONTINUE

      IF(MCMCflag.eq.1)then
       DO 106 J=1,NN
        DUSTH(J)=MCMCheight(J)
106    CONTINUE
      ENDIF

      IF(J.LT.2)THEN
        WRITE(*,*)' RDDMOD.f :: error reading dust-profile.'
        GOTO 119
      ENDIF
C Now sorting ...
      DO 114 K=1,NN
        DO 115 I=1,NN-1
          IF(ABS(DUSTH(I)-DUSTH(I+1)).LT.0.01)THEN
            WRITE(*,*)' RDDMOD.f :: identical heights found.'
            WRITE(*,*)K,I,DUSTH(I),DUSTH(I+1)
            GOTO 119
          ENDIF
          IF(DUSTH(I).GT.DUSTH(I+1))THEN
            WRITE(*,*)' RDDMOD.f :: reordering ...'
            WRITE(*,*)I,DUSTH(I),DUSTH(I+1)
            TD1 = DUSTH(I+1)
            DUSTH(I+1) = DUSTH(I)
            DUSTH(I) = TD1
            DO 233 L=1,NCONT
              TD1 = DUST(L,I+1)
              DUST(L,I+1) = DUST(L,I)
              DUST(L,I) = TD1
233         CONTINUE
          ENDIF
115     CONTINUE
114   CONTINUE

C Now interpolating the input array to find values at profile heights
      DO 244 K=1,NCONT
        DO 243 I=1,NN
          TDUST(I)=DUST(K,I)
243     CONTINUE
        DO 107 I=1,NPRO
          CALL VERINT(DUSTH,TDUST,NN,DUST(K,I),H(I))
107     CONTINUE
244   CONTINUE

119   CLOSE(UNIT=ilun)

      RETURN

      END
