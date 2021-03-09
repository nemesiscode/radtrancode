      SUBROUTINE NPARACON_ALL(vwave,p,t,ngas,idgas,isogas,amount,pp,
     1 fpara, xlen,absorb,iabsorb,dabsorb,idump)
C***********************************************************************
C_TITL:	NPARACON
C
C_DESC:	Computes a polynomial approximation to any known continuum
C	spectra for a variety of gas pairs over a defined wavenumber
C	region.
C
C_ARGS:	Input variables:
C	vwave		REAL		Wavenumber.
C	p		REAL		Total pressure [atm].
C	t		REAL		Temperature [K].
C	ngas		INTEGER		Number of gases.
C	idgas(ngas)	INTEGER		Array of gas identifiers.
C	isogas(ngas)	INTEGER		Array of gas isotope identifiers.
C	amount(ngas)	REAL		Amount of each gas 
C					[molecules/cm2].
C	pp(ngas)	REAL		Partial pressure of each gas
C	fpara		REAL		Para-H2 fraction.
C	xlen		REAL		path length in km
C	idump		INTEGER		If =1 then print diagnostic print
C					statements.
C
C	Output Variables:
C 	absorb          REAL		Calculated optical depth.
C	iabsorb(5)	INTEGER		Flag set to gas number 
C					corresponding to calculated gas
C					gradients listed below is active.
C	  DABSORB(7)	  Real	Rate of change of optical depth with:
C			  	1: H2 vmr
C			  	2: He vmr
C			        3: N2 vmr
C			        4: CH4 vmr
C				5: CO2 vmr
C			        6: Temperature
C			        7: para-H2 fraction
C
C_CALL:	fpread_all	Reads in the collision-induced absorption
C			coefficients of molecular hydrogen and helium from
C			0 to 1500 cm-1, for a range of different para-H2
C			fractions, as well as five other gases
C
C_HIST:	26nov86	SBC	ORIGINAL VERSION (HYDCON.F).
C	3feb88	SBC	Modified to return only one polynomial rather than
C			for all layers in GENLBL.
C	10feb97	CAN	This program adapted from HYDCON.F (H2-H2 and
C			H2-He continuum code) and its subroutine
C			OD_HYDTAB.F to create CIACON.F, which calculates
C			general gas pairs (mainly for Titan work).  
C	27jul01	PGJI	This version does gradients too! 
C	24oct03 PGJI	DH2DFP and DHEDFP (and hence DABSORB(6))
C			reinstated and debugged. We need these to compute
C			rate of change of H2-He CIA with para-H2 fraction
C	21oct05 PGJI	DABSORB array increased to length 7 to incorporate
C			CO2 CIA.
c     Fletcher (31/01/11) Added parameter DF to subroutine fpread.f
c	Fletcher (22/06/17) Updated to read in 5 additional gases
C***********************************************************************

        IMPLICIT NONE
        INCLUDE '../includes/constdef.f'
        INCLUDE '../includes/ciacom.f'
	INTEGER i,j,ngas
        INTEGER idgas(ngas),isogas(ngas),idump,iabsorb(5)
	REAL vwave,amount(ngas),p,t,sum,absorb,dabsorb(7)
        REAL xlen,height,height1,x1,f
        REAL PP(NGAS),tau,xp1
        REAL deltaf
        INTEGER KPARA

        REAL p0,t0,qh2,qhe,qh2s
C QH2: Molecular hydrogen mixing ratio.
C QHE: Helium mixing ratio.
        PARAMETER (p0 = 1.0,t0 = 273.15)
C Standard temperature and pressue (STP) equals 273.15 K, 1 atm.

        REAL MODBOLTZA,amag1,totam
        PARAMETER (MODBOLTZA = 10.*KBOLTZMANN/1.013)
C MODBOLTZ = KBOLTZ/1.013 (where KBOLTZ = 1.381E-23) and multiplied
C by 10. Then AMOUNT*MODBOLTZ*T(K)/P(ATM) gives path length in cm.
C AMAGAT = Number density at STP

        REAL FPARA,H2ALPHA(NPARA),HEALPHA(NPARA),AH2H2,AH2HE
        REAL DAH2DT(NPARA),DAHEDT(NPARA),DH2DT,DHEDT,DH2DFP,DHEDFP	
	real alpha(5),dalphadT(5),co2cia, aco2,n2n2cia
        real n2h2cia,an2n2,an2h2
        real qch4, qn2, qco2
        integer idiag,iquiet
        common/diagnostic/idiag,iquiet
C********************************* CODE ********************************

C       the mixing ratios
        qh2=0.
        qhe=0.
        qn2=0.
        qch4=0.
        qco2=0.
	
	
        DO i = 1, 5
	  IABSORB(i) = 0
	  DABSORB(i) = 0.0
	ENDDO
	DABSORB(6) = 0.0
	DABSORB(7) = 0.0

C-----------------------------------------------------------------------
C
C	First of all, determine which gases are present in the profile
C	and get their mixing ratios; the actual identifying numbers are
C	hard-wired (see the Radtran guide).
C
C-----------------------------------------------------------------------


        DO 30 i = 1, ngas
          IF ((idgas(i).EQ.39).AND.((isogas(i).EQ.0)
     &    .OR.(isogas(i).EQ.1))) THEN
            qh2 = pp(i)/p
            iabsorb(1) = i
	  ENDIF
          IF (idgas(i).EQ.40) THEN
            qhe = pp(i)/p
            iabsorb(2) = i
	  ENDIF
	  if (IDGAS(i) .eq. 22) then
		qn2 = pp(i) /p
		iabsorb(3)=i
	  endif
          if ( (IDGAS(i) .eq. 6) .and. 
     &     ( (ISOGAS(i) .eq. 0) .or. (ISOGAS(i) .eq. 1) ) ) then
          	qch4 = (pp(i)/p)
		iabsorb(4)=i
          endif
          if (IDGAS(i) .eq. 2) then
             qco2 = pp(i) /p
             iabsorb(5)=i
          endif
30      CONTINUE

C	Diagnostic - save the qh2 value if setting to zero later...
C	qh2s=qh2

C-----------------------------------------------------------------------
C
C	Calculate the height of the atmospheric layer (can use any one gas
C	to calc.)
C
C       height       ...... height in cm
C       amount(ngas) ...... no. of molecules per cm2 for a gas
C
C-----------------------------------------------------------------------

        IF (idump.EQ.1) WRITE(*,*)' NPARACON_ALL.f :: xlen = ',xlen

        IF (xlen.EQ.0) THEN
          xp1 = 0.0
          DO i = 1, ngas
            IF (pp(i).GT.XP1) THEN
              height = amount(i)*MODBOLTZA*t/pp(i)
              IF (idump.EQ.1) WRITE(*,*)' height (cm) = ',height
              xp1 = pp(i)
            ENDIF
          ENDDO      

C Calculate total path (length*(amagat)^2)
          tau = ((p/p0)**2)*((t0/t)**2)*height
          IF (idump.EQ.1) WRITE(*,*)' NPARACON.f :: tau = ',tau
        ELSE
          height = xlen*1e5
          IF (idump.EQ.1) WRITE(*,*)' height (cm) = ',height
          xp1 = -1.0
          DO I = 1, ngas
            IF(idump.eq.1)WRITE(*,*)'i,pp(i),xp1',i,pp(i),xp1
            IF (pp(i).GT.xp1.AND.pp(i).GT.0.0) THEN
              height1 = amount(i)*MODBOLTZA*t/pp(i)
              totam = amount(i)*p/pp(i)
              IF (idump.EQ.1) WRITE(*,*)' totam (cm-2) = ',TOTAM
              IF (idump.EQ.1) WRITE(*,*)' height1 (cm) = ',height1
              xp1 = pp(i)
            ENDIF
          ENDDO
          amag1 = (totam/height)/amagat
          tau = height*amag1**2
          IF (idump.EQ.1) WRITE(*,*)' NPARACON_ALL.f :: tau = ',tau
        ENDIF

        IF (idump.EQ.1) WRITE(*,*)' NPARACON_ALL.f :: qh2,qhe,qch4 = ',
     &	qh2,qhe,qch4
        IF (idump.EQ.1) WRITE(*,*)' NPARACON_ALL.f :: i, vwave: ',
     &	I,vwave


C-----------------------------------------------------------------------

	CALL FPREAD_ALL(vwave,T,KPARA,H2ALPHA,HEALPHA,DAH2DT,DAHEDT,
     1   deltaf,alpha,dalphadT)


c        CALL FPREAD(vwave,T,KPARA,H2ALPHA,HEALPHA,DAH2DT,DAHEDT,
c     1   deltaf)


C-----------------------------------------------------------------------
C	Now compute the total absorption.
C 	Start with H2-H2 and H2-He with the different para-H2 fractions

        IF (fpara.LT.0.0) THEN
          AH2H2 = H2ALPHA(KPARA)
          AH2HE = HEALPHA(KPARA)
          DH2DT = DAH2DT(KPARA)
          DHEDT = DAHEDT(KPARA)
          DH2DFP = 0.0
          DHEDFP = 0.0
        ELSE
          x1 = 1.0 + (fpara - 0.25)/deltaf 
          j = INT(x1)
          f = x1 - j
C          print*,'NPARACON: x1,deltaf,fpara,j,f: ', x1,deltaf,fpara,j,f
          IF (j.LT.1) THEN
            j = 1
            f = 0.0
          ELSE IF (j.GT.KPARA-2) THEN
            j = KPARA-2
            f = 1.0
          ENDIF

          AH2H2 = (1.0 - f)*H2ALPHA(j) + f*H2ALPHA(j+1)
          AH2HE = (1.0 - f)*HEALPHA(j) + f*HEALPHA(j+1)
          DH2DT = (1.0 - f)*DAH2DT(j) + f*DAH2DT(j+1)
          DHEDT = (1.0 - f)*DAHEDT(j) + f*DAHEDT(j+1)
          DH2DFP = (H2ALPHA(j+1) - H2ALPHA(j))/deltaf	!r.o.c.with paraH2
          DHEDFP = (HEALPHA(j+1) - HEALPHA(j))/deltaf
        ENDIF

C	Test of other gases - turn off H2 and He absorption
C	qh2=0.
C	qhe=0.
C	iabsorb(1)=0
C	iabsorb(2)=0

        SUM = AH2H2*qh2*qh2 + AH2HE*qh2*qhe
 	DABSORB(1) = (2*qh2*AH2H2 + qhe*AH2HE)
 	DABSORB(2) = qh2*AH2HE
 	DABSORB(6) = (DH2DT*qh2*qh2 + DHEDT*qh2*qhe)
 	DABSORB(7) = (DH2DFP*qh2*qh2 + DHEDFP*qh2*qhe)



C	Now add in contributions from H2-N2

	sum = sum + ( alpha(1) * qh2 * qn2 )
	DABSORB(1)=DABSORB(1) + qn2*alpha(1)
	DABSORB(3)=DABSORB(3) + qh2*alpha(1)
	DABSORB(6)=DABSORB(6) + qh2*qn2*dalphadT(1)
	if(IDUMP.EQ.1)then
	  print*, vwave, 'cm-1: h2n2   alpha (cm-1) is ', 
     &     (alpha(1) * qh2 *qn2)
	endif

C	Now add in contributions from N2-CH4

	sum = sum + ( alpha(2) * qn2 * qch4 )
	DABSORB(3)=DABSORB(3) + qch4*alpha(2)
	DABSORB(4)=DABSORB(4) + qn2*alpha(2)
	DABSORB(6)=DABSORB(6) + qn2*qch4*dalphadT(2)
	if(IDUMP.EQ.1)then
	print*, vwave, 'cm-1: n2ch4  alpha (cm-1) is ', 
     &     (alpha(2) * qn2 *qch4)
	endif

C	Now add in contributions from N2-N2

	sum = sum + ( alpha(3) * qn2 * qn2 )
	DABSORB(3)=DABSORB(3) + 2*qn2*alpha(3)
	DABSORB(6)=DABSORB(6) + qn2*qn2*dalphadT(3)
	if(IDUMP.EQ.1)then
	 print*, vwave, 'cm-1: n2n2   alpha (cm-1) is ', 
     &     (alpha(3) * qn2 *qn2)
	endif

C	Now add in contributions from CH4-CH4

	sum = sum + ( alpha(4) * qch4 * qch4 )
	DABSORB(4)=DABSORB(4) + 2*qch4*alpha(4)
	DABSORB(6)=DABSORB(6) + qch4*qch4*dalphadT(4)
	if(IDUMP.EQ.1)then
	 print*, vwave, 'cm-1: ch4ch4 alpha (cm-1) is ', 
     &     (alpha(4) * qch4 *qch4)
	endif

C	Now add in contributions from H2-CH4
C	Diagnostic, if all qh2 turned off except this one:
C	qh2=qh2s
	
	
	sum = sum + ( alpha(5) * qh2 * qch4 )
	DABSORB(1)=DABSORB(1) + qch4*alpha(5)
	DABSORB(4)=DABSORB(4) + qh2*alpha(5)
	DABSORB(6)=DABSORB(6) + qh2*qch4*dalphadT(5)
	if(idump.eq.1)then
	 print*, vwave, 'cm-1: ch4h2 alpha (cm-1) is ', 
     &     (alpha(5) * qh2 * qch4)
	endif

C	Look up CO2-CO2 CIA coefficients
	aco2 = co2cia(vwave)
	sum = sum + aco2 * qco2 * qco2
	DABSORB(5)=DABSORB(5) + 2*qco2*aco2

C	Look up N2-N2 NIR CIA coefficients
	an2n2 = n2n2cia(vwave)
	sum = sum + an2n2 * qn2 * qn2
	DABSORB(3)=DABSORB(3) + 2*qn2*an2n2

C	Look up N2H2 NIR CIA coefficients
	an2h2 = n2h2cia(vwave)
	sum = sum + an2h2 * qn2 * qh2
	DABSORB(3)=DABSORB(3) + qh2*an2h2
	DABSORB(1)=DABSORB(1) + qn2*an2h2


 	ABSORB = sum*tau

 	DO I=1,7
 	 DABSORB(I)=DABSORB(I)*tau
 	ENDDO



        IF (idump.EQ.1) THEN
          WRITE(*,*)vwave, 'cm-1: total (alpha*tau) is ',ABSORB
        ENDIF

	RETURN

	END
