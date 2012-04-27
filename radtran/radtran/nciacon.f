      SUBROUTINE NCIACON(V0,P,T,INORMAL,NGAS,IDGAS,
     &     ISOGAS, AMOUNT, PP, XLEN, ABSORB, IABS, DABSORB,
     &     IDUMP)
C----------------------------------------------------------------------------
C_TITLE:  CIACON: to compute gaseous continuum spectra from
C         a variety of gas pairs.
C
C Input variables:
C         V0:REAL         wavenumber
C         P:REAL          total pressure in atm
C         T:REAL          temp in K
C         INORMAL:INT     flag for ortho:para ratio (0=equilib. =1:1)
C                                                   (1=normal   =3:1)
C         NGAS:INT        no. of gases
C         IDGAS:INT(NGAS)  array of gas identifiers
C         ISOGAS:INT(NGAS) array of gas isotope identifiers
C         AMOUNT:REAL(NGAS)  amount of each gas (no. molecules/cm2)
C	  PP:REAL(NGAS)	  partial pressure of each gas
C	  XLEN:REAL 	  path length in km
C	  IDUMP:INT	  Set to  1 to print diagnostics
C
C Output Variables
C	  ABSORB	  Real	Calculated optical depth
C	  IABS(5)	  Int	Flag set to gas number corresponding to 
C				to calculated gas gradients listed below
C				gas gradients listed below is active
C	  DABSORB(7)	  Real	Rate of change of optical depth with:
C			  	1: H2 vmr
C			  	2: He vmr
C			        3: N2 vmr
C			        4: CH4 vmr
C				5: CO2 vmr
C			        6: Temperature
C			        7: para-H2 fraction
C_KEYS:   SUBR,ATMO,SPEC,VMS
C
C_DESCR:  computes a polynomial approximation to any known continuum spectra
C         for a particular gas over a defined wavenumber region.
C
C_FILES:  none
C
C_CALLS   none
C
C_BUGS:
C
C_HIST:   26nov86 SBC ORIGINAL VERSION (HYDCON.F)
C          3feb88 SBC modified to return only one polynomial rather than
C                 for all layers in GENLBL
C         10feb97 This program adapted from HYDCON.F (H2-H2 and H2-He 
C                 continuum code) and its subroutine OD_HYDTAB.F
C                 to create CIACON.F, which calculates general gas pairs
C                 (mainly for Titan work) (C.Nixon)
C	  27jul01 This version does gradients too!
C
C----------------------------------------------------------------------------
        INCLUDE '../includes/constdef.f'
        INCLUDE '../includes/ciacom.f'
	INTEGER ID,ISO,I,J, INORMAL, NGAS
        INTEGER IDGAS(NGAS), ISOGAS(NGAS), IDUMP, IABS(5)
	REAL V0,AMOUNT(NGAS),p,T,ABSORB, sum, DABSORB(7)
        real xlen,height,height1
        REAL XMIN,XMAX,MODBOLTZA,PP(NGAS),tau,z
        real P0,T0,amag1
        parameter (P0=1.,T0=273.15)
        PARAMETER (MODBOLTZA=10.*KBOLTZMANN/1.013)
C	MODBOLTZA = KBOLTZ/1.013 (where KBOLTZ = 1.381E-23) and multiplied
C 	by 10. Then AMOUNT*MODBOLTZ*T(K)/P(ATM) gives path length in cm.
C	AMAGAT = Number density at STP (1 atm, 273.15 K) 
        real alpha(NUMPAIRS),dalphadT(NUMPAIRS),co2cia, aco2,n2n2cia
        real n2h2cia,an2n2,an2h2
        real qhe, qh2, qch4, qn2, qco2

C	In IABSORB(5)	1=H2
C			2=He
C			3=N2
C			4=CH4
C		 	5=CO2
C	DABSORB(1-5) is as above. DABSORB(6) is dtaudtemp, DABSORB(7) is 
C	dtaudfpara (only set in nparacon.f)
C

C       the mixing ratios
        qh2=0.
        qhe=0.
        qn2=0.
        qch4=0.
        qco2=0.

        DO I=1,5
         IABS(I)=0
	 DABSORB(I)=0.0
        ENDDO
        DABSORB(6)=0.0
        DABSORB(7)=0.0

C        print*,'nciacon'
C        print*,V0,P,T,INORMAL,NGAS
C        print*,(IDGAS(i),i=1,ngas)
C        print*,(ISOGAS(i),i=1,ngas)
C        print*,(AMOUNT(i),i=1,ngas)
C        print*,(PP(i),i=1,ngas)
C        print*,IDUMP
C        print*,(IABS(i),i=1,5)
C        print*,(DABSORB(i),i=1,7)
 

C       first of all, figure out which gases are present in the profile
C       and get their mixing ratios
C       the actual identifying numbers are hard-wired (see Radtran guide)
        do 30 i=1, NGAS


           if ( (IDGAS(i) .eq. 39) .and.
     &     ( (ISOGAS(i). eq. 0) .or. (ISOGAS(i) .eq. 1) ) ) then
          	qh2 = pp(i) /p
		iabs(1)=i
	   endif
           if (IDGAS(i) .eq. 40) then
		qhe = pp(i) /p
		iabs(2)=i
	   endif
           if (IDGAS(i) .eq. 22) then
		qn2 = pp(i) /p
		iabs(3)=i
	   endif
           if ( (IDGAS(i) .eq. 6) .and. 
     &     ( (ISOGAS(i) .eq. 0) .or. (ISOGAS(i) .eq. 1) ) ) then
          	qch4 = (pp(i)/p)
		iabs(4)=i
           endif
           if (IDGAS(i) .eq. 2) then
             qco2 = pp(i) /p
             iabs(5)=i
           endif
 30     continue


        
C       calculate the height of the atmospheric layer
C       (can use any one gas to calc.)
C
C       height       ...... height in cm
C       amount(ngas) ...... no. of molecules per cm2 for a gas

        IF(IDUMP.EQ.1)print*,'NCIACON. xlen = ',xlen

        if(xlen.eq.0)then
         XP1 = 0.0
         do I=1,NGAS
          if(pp(i).gt.XP1)then
            height = AMOUNT(i)*MODBOLTZA*T/pp(i)
            if(IDUMP.EQ.1)print*,'height (cm) = ',height
            XP1=pp(i)
          endif
         enddo      

C        Calculate total path length*(amagat)^2
         tau=((P/P0)**2)*((T0/T)**2)*height
         if(IDUMP.eq.1)print*,'tau = ',tau

        else

         height = xlen*1e5
         if(IDUMP.EQ.1)print*,'height (cm) = ',height
         XP1 = -1.0
         do I=1,NGAS
          if(pp(i).gt.XP1.and.pp(i).gt.0.0)then
            height1 = AMOUNT(i)*MODBOLTZA*T/pp(i)
            TOTAM = AMOUNT(i)*P/pp(i)
            if(IDUMP.EQ.1)print*,'totam (cm-2) = ',TOTAM
            if(IDUMP.EQ.1)print*,'height1 (cm) = ',height1
            XP1=pp(i)
          endif
         enddo      
         amag1 = (TOTAM/height)/AMAGAT
         tau = height*amag1**2
         if(IDUMP.eq.1)print*,'tau = ',tau
        endif

        if(IDUMP.EQ.1)print*,'nciacon. qh2,qhe,qch4 = ',qh2,qhe,qch4


        if(IDUMP.EQ.1)print*,'V0 : ',V0
C       this will hold the total opacity
        sum = 0.

C       goes away and returns the abs. co-eff for all gases in the table
C       at a given t and wavenumber
        call ciaread(V0,T,alpha,dalphadT)

C       to use the `alpha's we need to know the ordering.
C       nine gas pairs at present:
C
C 1............H2-H2 (ortho:para = 1:1 `equilibrium')
C 2............H2-He                        "
C 3............H2-H2 (ortho:para = 3:1 `normal')
C 4............H2-He                        "
C 5............H2-N2
C 6............N2-CH4
C 7............N2-N2
C 8............CH4-CH4
C 9............H2-CH4 
C 10...........CO2-CO2


C       equilibrium hydrogen
        if (INORMAL .eq. 0) then
           sum = sum + ( alpha(1) * qh2 * qh2 ) +
     &		( alpha(2) * qh2 * qhe )
	   DABSORB(1)=DABSORB(1) + 2*qh2*alpha(1) + qHe*alpha(2)
	   DABSORB(2)=DABSORB(2) + qh2*alpha(2)
	   DABSORB(6)=DABSORB(6)+qh2*qh2*dalphadT(1) + 
     &		qh2*qhe*dalphadT(2)

           if(IDUMP.eq.1)then
            print*, V0, 'cm-1: h2h2   alpha (cm-1) is ', 
     &     (alpha(1) * qh2 *qh2)
            print*, V0, 'cm-1: h2he   alpha (cm-1) is ', 
     &     (alpha(2) * qh2 *qhe)
           endif
        else
C `normal' hydrogen
           sum = sum + ( alpha(3) * qh2 * qh2) +
     &			( alpha(4) * qh2 * qhe )
	   DABSORB(1)=DABSORB(1) + 2*qh2*alpha(3) + qhe*alpha(4)
	   DABSORB(2)=DABSORB(2) + qh2*alpha(4)
	   DABSORB(6)=DABSORB(6)+qh2*qh2*dalphadT(3) + 
     &		qh2*qhe*dalphadT(4)
           if(IDUMP.EQ.1)then
            print*, V0, 'cm-1: h2h2   alpha (cm-1) is ', 
     &     (alpha(3) * qh2 *qh2)
            print*, V0, 'cm-1: h2he   alpha (cm-1) is ', 
     &     (alpha(4) * qh2 *qhe)
           endif
        end if

        sum = sum + ( alpha(5) * qh2 * qn2 )
        DABSORB(1)=DABSORB(1) + qn2*alpha(5)
	DABSORB(3)=DABSORB(3) + qh2*alpha(5)
	DABSORB(6)=DABSORB(6) + qh2*qn2*dalphadT(5)
        if(IDUMP.EQ.1)then
          print*, V0, 'cm-1: h2n2   alpha (cm-1) is ', 
     &     (alpha(5) * qh2 *qn2)
        endif

        sum = sum + ( alpha(6) * qn2 * qch4 )
  	DABSORB(3)=DABSORB(3) + qch4*alpha(6)
	DABSORB(4)=DABSORB(4) + qn2*alpha(6)
	DABSORB(6)=DABSORB(6) + qn2*qch4*dalphadT(6)
        if(IDUMP.EQ.1)then
        print*, V0, 'cm-1: n2ch4  alpha (cm-1) is ', 
     &     (alpha(6) * qn2 *qch4)
        endif

        sum = sum + ( alpha(7) * qn2 * qn2 )
  	DABSORB(3)=DABSORB(3) + 2*qn2*alpha(7)
	DABSORB(6)=DABSORB(6) + qn2*qn2*dalphadT(7)
        if(IDUMP.EQ.1)then
         print*, V0, 'cm-1: n2n2   alpha (cm-1) is ', 
     &     (alpha(7) * qn2 *qn2)
        endif

        sum = sum + ( alpha(8) * qch4 * qch4 )
  	DABSORB(4)=DABSORB(4) + 2*qch4*alpha(8)
	DABSORB(6)=DABSORB(6) + qch4*qch4*dalphadT(8)
        if(IDUMP.EQ.1)then
         print*, V0, 'cm-1: ch4ch4 alpha (cm-1) is ', 
     &     (alpha(8) * qch4 *qch4)
        endif

        sum = sum + ( alpha(9) * qh2 * qch4 )
  	DABSORB(1)=DABSORB(1) + qch4*alpha(9)
  	DABSORB(4)=DABSORB(4) + qh2*alpha(9)
	DABSORB(6)=DABSORB(6) + qh2*qch4*dalphadT(9)
        if(idump.eq.1)then
         print*, V0, 'cm-1: ch4h2 alpha (cm-1) is ', 
     &     (alpha(9) * qh2 * qch4)
        endif

C       Look up CO2-CO2 CIA coefficients
        aco2 = co2cia(V0)
        sum = sum + aco2 * qco2 * qco2
        DABSORB(5)=DABSORB(5) + 2*qco2*aco2

C       Look up N2-N2 NIR CIA coefficients
        an2n2 = n2n2cia(V0)
        sum = sum + an2n2 * qn2 * qn2
        DABSORB(3)=DABSORB(3) + 2*qn2*an2n2

C       Look up N2H2 NIR CIA coefficients
        an2h2 = n2h2cia(V0)
        sum = sum + an2h2 * qn2 * qh2
        DABSORB(3)=DABSORB(3) + qh2*an2h2
        DABSORB(1)=DABSORB(1) + qn2*an2h2


        ABSORB = sum*tau

        DO I=1,7
	  DABSORB(I)=DABSORB(I)*tau
        ENDDO

        if(IDUMP.eq.1)then
          print*, V0, 'cm-1: total (alpha*tau) is ',ABSORB
        endif
  
	RETURN
	END

