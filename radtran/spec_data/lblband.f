      SUBROUTINE LBLBAND(POUT,WAVEN)
C     $Id: lblband.f,v 1.5 2011-06-17 15:53:01 irwin Exp $
C-------------------------------------------------------------------------
C_TITLE: LBLBAND: Calculates the Goody-Voigt or Malkmus-Voigt equivalent 
C        band parameters by statistical analysis of the linedata. 
C        Goody-Lorentz and Malkmus-Lorentz lines are easily extracted from
C        these. 
C
C_ARGS:  NGAS:INTEGER     number of gases to consider
C        IDGAS(MAXGAS):INTEGER   the LOCAL gas identifier for each gas
C                         The local gas identifier agrees with HITRAN id's
C                         as far as possible (Jan 1988) with extra id's
C                         for gases not in the HITRAN compilation. eg.
C                         those in the GEISA compilation
C        ISOGAS(MAXGAS):INTEGER   the local isotopic identifier, if zero all
C                         isotopes of the gas are included.
C                         Isotope id's also agree as far as possible with
C                         HITRAN id's. Similar (1,2,3...n) id's have been
C                         defined for additional gases.                    
C                         If zero then line strengths are used as tabulated
C                         (i.e. corrected for normal terrestrial distribution)
C                         If a specific isotope is selected then its strength
C                         is scaled to the pure isotope using the HITRAN
C                         data base values
C        VMIN:REAL        minimum wavenumber for output
C        DELV:REAL        wavenumber spacing for output
C        NPOINT:INTEGER   the number of output points
C                         i.e. each wavenumber, v=vmin+(i-1)*delv for
C                         i=1 to npoint
C
C_FILES  DBLUN - line data base files
C        none unless used by a particular OUTPUT model. See code for
C        descriptions of model dependent parameters and input
C
C_CALLS: 
C
C_BUGS:
C
C_HIST:  23nov93 PGJI Original version. Adapted from SBC's GENLBL 
C        23nov93 PGJI Stimulated emission term included
C 	 02jun94 PGJI Converted from GENLBL for band parameter applications
C	 02jun94 PGJI Calculates the statistical band parameters required by
C     		      Kim's formulation of the Goody-Voigt Model i.e.:
C		      Knu0,delta/AD,AL/AD,EL and SFB
C	   	      For Goody-Lorentz model AL is the same (AL/AD * AD/delta)
C		      For Malkmus-Lorentz model AL is 0.25*this value
C        30nov95 PJCS Converted AD0 to an array to stop the mass of the final
C                     gas being used for all the gases.
C------------------------------------------------------------------------------
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/lincom.f'
      INCLUDE '../includes/dbcom.f' 
C------------------------------------------------------------------------------
      REAL POUT(MAXOUT,MAXGAS,7),WAVEN(MAXOUT,2)
      REAL AD,V0,V,V1,DELWAV
      INTEGER GAS,NGAS1,IDGAS1(1),ISOGAS1(1)
      REAL WR,SR,SSFB,SEL,ABSCO,AL,SFB,EL
      REAL DELTA,KNU0,Y0,DELAD,PI
      PARAMETER (PI=3.1415927)
      DOUBLE PRECISION DX0,DX1
      INTEGER FSTLIN,LSTLIN,NLIN
      CHARACTER*1 ANS
      DIMENSION AD0(MAXGAS)

C------------------------------------------------------------------------------
C     code section
C----------
      ISOMAX=MAXISO
C     ISOMAX used in subroutine calls - maximum number of isotopes that
C     can be used (per gas)
      LINMAX=MAXLIN

      DO 16 I=1,NGAS
C     checking isotope included in model and setting mass for doppler
C     width calculation
      IF(ISOGAS(I).EQ.0)THEN
C       have to compute mass for normal combination if using all isotopes
        XMASS=0.
        DO 255 J=1,DBNISO(IDGAS(I))
        XMASS=XMASS+MASSNO(J,IDGAS(I))*RELABU(J,IDGAS(I))
255     CONTINUE
C       it is possible that the relative abundance for one of the isotopes is
C       wrong (eg set to 1 because unknown) so checking that the final
C       mass is within 20% of the mass of the main isotope
        IF(MASSNO(1,IDGAS(I)).GT.1.E-32.AND.
     1  ABS((XMASS-MASSNO(1,IDGAS(I)))/MASSNO(1,IDGAS(I)))
     1  .GT.0.2)THEN
         XMASS=MASSNO(1,IDGAS(I))
         WRITE(*,287)
287      FORMAT(' %warning - using main isotope mass number')
         END IF
      ELSE
        IF(ISOGAS(I).LE.DBNISO(IDGAS(I)))GOTO 418
C       isotopes are included in arrays in same order as new HITRAN 86
C       definitions 
        WRITE(*,417)ISOGAS(I),IDGAS(I)
417     FORMAT(' %model doesn"t include isotope ',I5,' for gas ',I2)
        STOP
418     CONTINUE
        XMASS=MASSNO(ISOGAS(I),IDGAS(I))
C
      END IF
      WRITE(*,420)IDGAS(I),ISOGAS(I),GASNAM(IDGAS(I)),XMASS
420   FORMAT(' %gas:',I2,' / ',I2,2X,1A6,' mass=',F7.2)
      IF(ISOGAS(I).EQ.0)THEN
        K=1
       ELSE
        K=ISOGAS(I)
      END IF
      AD0(I)=4.301E-7/SQRT(XMASS)
16    CONTINUE

C     First read in the linedata for the bin

      DO 5000 I=1,NPOINT

      V0=VMIN+DELV*(I-1)-0.5*FWHM
      V=V0 + 0.5*FWHM
      V1=V0+FWHM
      DELWAV=FWHM

C     Status line to help user keep track of progress
      WRITE(*,430) V
 430  FORMAT(' V = ',F8.2)
C     End of status line      

      FSTLIN=1
      DX0=DBLE(V0)
      DX1=DBLE(DELWAV)
      CALL LINESS(DX0,DX1,MAXLIN,VLIN,SLIN,ALIN,ELIN,IDLIN,SBLIN,
     & PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,FSTLIN,LSTLIN,NGAS,IDGAS,ISOGAS)

      print*,'Total lines read in bin = ',I,FSTLIN,LSTLIN,NLIN

      DO 52 GAS=1,NGAS
C        print*,gas
        AD=AD0(GAS)*V

        WR=0.
        SR=0.
        SSFB=0.
        SEL=0.

        DO 55 LINE=1,NLIN
         IF(IDLIN(LINE).EQ.GAS)THEN    
C         Extra section to stop very small strengths coming in. Here limit
C         to S > 1e-36
C          IF(SLIN(LINE).GE.1.0E-9)THEN
C           ABSCO=SLIN(LINE)*1E-27
C          ELSE
C           ABSCO = 0.
C          END IF
          ABSCO=SNGL(SLIN(LINE)*1E-27)
          AL = ALIN(LINE)-SBLIN(LINE)
          SFB = (ALIN(LINE)-SBLIN(LINE))/ALIN(LINE)
          EL = ELIN(LINE)
C         Test to discard lines with negative LSE

          print*,VLIN(LINE),ABSCO,AL,EL,SFB
          IF(EL.LT.0.0)ABSCO=0.0
          IF(ABSCO.GT.0)THEN
           WR=WR+ABSCO
           SR=SR+SQRT(ABSCO*AL)
           SSFB=SSFB+SFB*ABSCO
           SEL=SEL+EL*ABSCO
          ENDIF
          print*,WR,SR,SSFB,SEL
          read(5,1)ans
1         format(a)
         ENDIF
55      CONTINUE
 
        IF(WR.GT.0.)THEN
         DELTA = DELWAV/REAL(NLIN)
         SFB=SSFB/WR
         EL=SEL/WR
         WR=WR/DELWAV
         SR=2.*SR/DELWAV
         KNU0=WR
         AL=DELTA*SR*SR/(PI*KNU0)
         Y0 = AL/AD
         DELAD = DELTA / AD
        ELSE
         KNU0=0.
         DELAD=1.
         Y0=1.
         EL=1.
         SFB=1. 
        ENDIF

C       Begining of debugging lines
        PRINT*,'  V    GAS    Kv(T0)      delta/AD0         y0
     &            El            SFB'
        WRITE(*,3267) V,GAS,KNU0,DELAD,Y0,EL,SFB
3267    FORMAT(F7.1,I3,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
C       End of debugging lines

        read(5,1)ans

        POUT(I,GAS,1) = KNU0
        POUT(I,GAS,2) = DELAD
        POUT(I,GAS,3) = Y0
        POUT(I,GAS,4) = EL
        POUT(I,GAS,5) = SFB

        WAVEN(I,1)=V
        WAVEN(I,2)=DELWAV

52    CONTINUE
      
5000  CONTINUE

      RETURN
      END
C


