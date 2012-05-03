      subroutine clouds(npro,h,P,T,deepH2O,deepNH3,deepH2S,
     &                   ssatH2O, ssatNH3,
     &                   cbaseH2O,cbaseNH4SH,cbaseNH3,
     &                   shH2O,shNH4SH,shNH3,vmrH2O,vmrNH3,vmrH2S)
C     $Id: clouds.f,v 1.4 2011-01-20 14:30:32 teanby Exp $
C     *******************************************************************
C     This is a fortran implementation of a routine based on the cloud 
C     model written by Philip Smith based on Lewis's model.
C     *******************************************************************
C     Input variables:
C	npro	integer	number of vertical ordinates
C	h(npro)	real	height ordinates (km)
C	P(npro)	real	pressure ordinates (atm)
C	T(npro)	real	temperature ordinates (K)
C	deepH2O	real	Deep H2O volume mixing ratio
C	deepNH3	real	Deep NH3 "      "      "
C	deepH2S	real	Deep H2S "      "      "
C       ssatH2O real    Level of supersaturation of H2O _above_ H2O cloud base. 
C                         1.1=10% supersaturation
C       ssatNH3 real    Level of supersaturation of NH3 _above_ NH3 cloud base. 
C                         1.1=10% supersaturation
C       NB: it is not clear what might be meant by supersaturtion for the NH4SH cloud,
C           so it is currently ignored.
C
C     Output variables
C	cbaseH2O	real	Cloud base height (H2O) in km
C	cbaseNH4SH	real	Cloud base height (NH4SH) in km
C	cbaseNH3	real	Cloud base height (NH3) in km
C       shH2O	        real	Scale height (H2O) in km
C       shNH4SH	        real	Scale height (NH4SH) in km
C	shNH3	        real	Scale height (NH3) in km
C	vmrH2O(npro)	real	H2O volume mixing ratio profile
C	vmrNH3(npro)	real	NH3 volume mixing ratio profile
C	vmrH2S(npro)	real	N2S volume mixing ratio profile
C
C     *******************************************************************

      INTEGER  NPRO,MAXPRO,NUMGAS,NUMCON
      INTEGER RESLEVEL,RESINC
      PARAMETER (MAXPRO=200)
      PARAMETER (NUMGAS=3)   !This is the number of gaseous species
      PARAMETER (NUMCON=4)   !This is the number of condensible species
      PARAMETER (RESINC=20 ) !The height resolution is increased by 2^resinc for the cloud base

      REAL   h(NPRO), P(NPRO), T(NPRO)
      REAL   deepH2O, deepNH3, deepH2S
      REAL   ssatH2O, ssatNH3
      REAL   cbaseH2O, cbaseNH4SH, cbaseNH3
      REAL   shH2O, shNH4SH, shNH3
      REAL   vmrH2O(NPRO),vmrNH3(NPRO),vmrH2S(NPRO)
      REAL   vmr(NUMGAS,MAXPRO)  !vmr contains H2O, NH3, H2S in that order
      REAL   ssat(NUMGAS)        !super saturation variables in array form

C     See below for definition of DHlg, DHsl, DHsg, TT, TP, R

      DOUBLE PRECISION DHlg(NUMCON),DHsl(NUMCON),DHsg(NUMCON)
      DOUBLE PRECISION TT(NUMCON),TP(NUMCON)
      DOUBLE PRECISION R

C     'cloud_lay' = the number of layers in which condensation
C        has taken place for each condensate.
C     'sh' contains the three scale heights
C     'cbase' contains the three cloud bases
C     'cloud_den' contains the cloud density for each condensible
C        at each height in units of vmr/km. To convert this to units
C	 that mean something, ie mols/m^3
C        one must take into consideration the shape of the column 
C        within which the rising air is constrained, and the size 
C        of the air parcel that is rising. In this model it is assumed
C        that the column is a 'pipe' of constant x-section, which means
C        that one should also divide by the x-sectional area (unknown!).
C        However we are saved because we only wish to calculate the 
C        scale-height where any constant factor with height in the 
C        forumla for cloud_den cancels out, ie the x-sectional area. 
C
C        If one assumes that the volume at each level is given by the 
C        volume of one mole of gas, then the cloud density will be 
C        approximately proportional to the step size which is obviously
C        undesirable.
C 
C        The moral of all this is think carefully before using.

      INTEGER  cloud_lay(NUMCON)
      REAL   sh(NUMCON),cbase(NUMCON),cloud_den(NUMCON,MAXPRO)

C     set up dummy variables
      INTEGER  level,con,gas ! height level, condensate, and gas indexes
      REAL   c_vp     ! condensate vapour pressure
      REAL   vmrBC    !vmr Below Cloud
      REAL   K         ! Equilibrium constant for NH4SH <==> NH3 + H2S
      REAL   pnh3,ph2s ! partial pressure of NH3 & H2S before NH4SH precipitation.
      REAL   x         ! number of moles of NH4SH precipitated.
      INTEGER  cbase_index(numgas) ! the index of the cloud base 'layer' for the simple gas ppt's

C     -------- Set up constants ---------------
C     DHlg = enthalpy going from liquid to gas (J/mol)
C     DHsl = enthalpy going from solid to liquid (J/mol)
C     DHsg = enthalpy going from solid to gas (J/mol) [=DHsl+DHlg]
C     TT   = Triple point temperature (K)
C     TP   = Triple point pressure (bar)
C     R    = Gas constant (J/mol/K)

C     Data for the gases at the triple point (American IOP handbook &
C     International Critical Tables). Th DHsl for NH4SH is unknown, and the
C     value quoted for DHlg is actually the known value for DHsg.
C                     H20        NH3        H2S       NH4SH
      DATA DHlg / 45.0491D+3, 25.359D+3, 19.54D+3,  46.025D+3 /
      DATA DHsl /  6.0082D+3,  5.652D+3,  2.377D+3,       0.   /
      DATA TT   /     273.16,    195.40,   187.61,     391.   /
      DATA TP   /   6.026D-3,  5.970D-2,    0.229,     51.96  /
      DATA R    / 8.31434 /

      DO con = 1,numcon
         DHsg(con) = DHlg(con) + DHsl(con)
      ENDDO

C     -------- End of constants ---------------

C     -------- Start of debugging lines ---------
C      PRINT*,NPRO
C      PRINT*,H(NPRO), P(NPRO), T(NPRO)
C      PRINT*,DEEPH2O, DEEPNH3, DEEPH2S
C      PRINT*,ssatH2O, ssatNH3
C      PRINT*,  CBASEH2O, CBASENH3, CBASENH4SH
C      PRINT*,  SHH2O, SHNH3, SHNH4SH
C      PRINT*,VMRH2O(NPRO),VMRNH3(NPRO),VMRH2S(NPRO)
C     -------- End of debugging lines ---------

C     -------- Set up arrays at bottom of profile -------
C     We assume that all of the H2S reacts with NH3 and that
C     there is more NH3 than H2S.

      vmr(1,1)=deepH2O
      vmr(2,1)=deepNH3
      vmr(3,1)=deepH2S

C     Convert ssat input variables into an array form
      ssat(1) = ssatH2O
      ssat(2) = ssatNH3
      ssat(3) = 1.        !ssat H2S is not currently an input parameter                          

C     Set up array for the indicies of the heights of the cloud bases
C     The indicies are set to npro+1, ie off the top of the profile so 
C     that if that gas doesn't crystalize in its pure form it won't have
C     its index recorded, and therefore no supersaturation effect will
C     be applied.
      DO gas=1,numgas
         cbase_index(gas)=npro+1
      ENDDO

C     ========== Start of Main program ============

      DO level=2,NPRO
         DO con=1,numcon  
C     'con' is the number of the condensible, 1=H2O, 2=NH3, 3=H2S, 4=NH4SH
C     We assume here that we are condensing into the solid phase.
C     This may well be wrong for H2O, but shouldn't be too far out,
C     and is a greatly simplifying assumption.

C ----- H2O, NH3, H2S clouds ---------
            IF (con.LE.3) then     ! ie H2O, NH3, H2S
            c_vp = TP(con)*EXP(-DHsg(con)/R*(1/T(level)-1/TT(con)))
            IF (c_vp.LT.(vmr(con,level-1)*P(level))) THEN
               vmr(con,level)=c_vp/P(level)
               cloud_den(con,level)=(vmr(con,level-1)-vmr(con,level))/
     &(h(level)-h(level-1))        ! See notes on cloud_den above.
               IF (cloud_den(con,level).LT.0) THEN
                  PRINT*,"WARNING (clouds): Negative cloud density"
               ENDIF
               cloud_lay(con)=cloud_lay(con)+1
C     *** Start of debugging lines ***
CDEBUG  	     IF (con.EQ.2) THEN
CDEBUG  		print*
CDEBUG  		print*,"level = ",level,",   con = ",con
CDEBUG  		print*,"P   = ",P(level)
CDEBUG  		PRINT*,"vmr(old) = ",vmr(con,level-1)
CDEBUG  		print*,"c_vp = ",c_vp
CDEBUG  		print*,"vmr*P = ",vmr(con,level-1)*P(level)
CDEBUG  		PRINT*,"vmr(new) = ",vmr(con,level)
CDEBUG  	     ENDIF
C     *** End of debugging lines ***

C  An iterative method is used to determine the height of the cloud base
C  more accurately than the provided height spacing allows. This is done by 
C  sampling the mid point of the layer and then creating a 'new' layer by using
C  either the lower or upper half depending on whether the mid-point is in cloud
C  or not. This bisection method should reliably home in on the cloud layer. The
C  number or iterations is given by resinc. 
C
C  Linear interpolation in ln P and T is used.
               IF (cloud_lay(con).EQ.1) THEN
                  vmrBC=vmr(con,level-1)     !vmr Below Cloud
                  P1 = P(level-1)
                  T1 = T(level-1)
                  h1 = h(level-1)
                  P2 = P(level)
                  T2 = T(level)
                  h2 = h(level)
                  DO reslevel=1,resinc
                     Pmid=SQRT(P1*P2) !This is identical to the mid point in log space
                     Tmid=(T1+T2)/2
                     hmid=(h1+h2)/2
                     c_vp = TP(con)*
     &                  EXP(-DHsg(con)/R*(1/Tmid-1/TT(con)))
                     IF (c_vp.LT.(vmrBC*Pmid)) THEN
                        P2=Pmid
                        T2=Tmid
                        h2=hmid
                     ELSE
                        P1=Pmid
                        T1=Tmid
                        h1=hmid
                     ENDIF
                  ENDDO
                  cbase(con)=(h1+h2)/2  ! NB this is not hmid.
C Record level of cloud bottom to allow compensation for supersaturaion later
                  cbase_index(con) = level    
               ENDIF
            ELSE
               vmr(con,level)=vmr(con,level-1)
            ENDIF

C --- NH4SH cloud ------
C from Lewis, Icarus, 1969
C K=p'(NH3) * p'(H2S) at equilibrium
C log10(K) = 14.82-(4705/T)
C therefore:
C K=[p(NH3)-x*P][p(H2S)-x*P]
C where x is the number of moles of ppt
C NB: melting point of NH4SH is 118C=391K from CRC h/book of chem & phys.

C NB: we use level _not_ level-1 for initial vmr since the vmr for level is
C     already defined for the NH3 & H2S cloud routines above, and one could overwrite
C     the value of vmr(level) already calculated, or use the gas twice, although
C     this latter possibility is less important as the 3 main cloud layers are
C     well separated.

         ELSE IF (con.EQ.4) THEN
            K=10.**(14.82-(4705./T(level)))
            pnh3=vmr(2,level)*P(level)
            ph2s=vmr(3,level)*P(level)

C NB: we use the alternative form of the eqn for roots of a quadratic because
C it is more accurate. See "numerical recipes" for more info.
            x=2*(pnh3*ph2s-k)/P(level)/
     &         ((pnh3+ph2s)+SQRT((pnh3+ph2s)**2-4*pnh3*ph2s+4*K))

            IF (x.GT.0) THEN

               vmr(2,level)=vmr(2,level)-x
               vmr(3,level)=vmr(3,level)-x

               cloud_den(con,level)=(x)/(h(level)-h(level-1)) ! See notes on cloud_den above.
               IF (cloud_den(con,level).LT.0) THEN
                  PRINT*,"WARNING (clouds): Negative cloud density"
               ENDIF
               cloud_lay(con)=cloud_lay(con)+1

C     *** Start of debugging lines ***
C            print*
C            print*,"level = ",level,",   con = ",con
C            print*,"P   = ",P(level)
C            Print*,"K          = ",K
C            Print*,"[NH3][H2S] = ",vmr(2,level)*vmr(3,level)*P(level)**2
C            Print*,"x = ",x
C     *** End of debugging lines ***

C  An iterative method is used to determine the height of the cloud base
C  more accurately than the provided height spacing allows. This is done by 
C  sampling the mid point of the layer and then creating a 'new' layer by using
C  either the lower or upper half depending on whether the mid-point is in cloud
C  or not. This bisection method should reliably home in on the cloud layer. The
C  number or iterations is given by resinc. 
C
C  Linear interpolation in ln P and T is used.
               IF (cloud_lay(con).EQ.1) THEN
                  P1 = P(level-1)
                  T1 = T(level-1)
                  h1 = h(level-1)
                  P2 = P(level)
                  T2 = T(level)
                  h2 = h(level)
                  DO reslevel=1,resinc
                     Pmid=SQRT(P1*P2) !This is identical to the mid point in log space
                     Tmid=(T1+T2)/2
                     hmid=(h1+h2)/2
C   Calculate x at mid point:
            K=10**(14.82-(4705./Tmid))
            x=((pnh3+ph2s)-SQRT((pnh3-ph2s)**2+K))/(2*Pmid)
C   If x>0 then we would get a ppt at the mid point.
                     IF (x.GT.0) THEN   
                        P2=Pmid
                        T2=Tmid
                        h2=hmid
                     ELSE
                        P1=Pmid
                        T1=Tmid
                        h1=hmid
                     ENDIF
                  ENDDO
                  cbase(con)=(h1+h2)/2  ! NB this is not hmid.
               ENDIF
            ENDIF
         ENDIF

C     Calculate scale height of cloud:

         IF (cloud_lay(con).EQ.3) THEN
            sh(con)=(h(level)-h(level-1))/
     &           LOG(cloud_den(con,level-1)/cloud_den(con,level))
C     *** Start of debugging lines ***
C       PRINT*,"cloud_density (level 1, gas=",con,") = ",
C     &cloud_den(con,level-2)
C       PRINT*,"cloud_density (level 2, gas=",con,") = ",
C     &cloud_den(con,level-1)
C       PRINT*,"cloud_density (level 3, gas=",con,") = ",
C     &cloud_den(con,level)
C      PRINT*,"Scale Height for (gas=",con,") = ",sh(con)
C     *** End of debugging lines ***
         ENDIF
         ENDDO
      ENDDO


C     --------- Supersaturate gasses ---------
C Because of limitations in the model it may be desirable to 
C allow gasses to be super(under)saturated from cloud bases upwards
C without affecting the main cloud results. This is handled here by
C multiplying the vmr of each gas by the level of supersaturation (ssat)
C where ssat=1.1 means 10% supersaturated. 
C
C It is not clear what a supersaturation for NH4SH means, and so it is
C simply ignored. For each gas the supersaturation multiplier is only 
C applied to levels above the base of the ppt formed by precipitating 
C that gas only, ie the H2O, NH3, H2S clouds.
C
C This method with dealing with supersaturation is unphysical, so it may 
C produce bizare results if one isn't careful.

      DO gas=1,numgas
         DO level=cbase_index(gas),npro
            vmr(gas,level)=ssat(gas)*vmr(gas,level)
         ENDDO
      ENDDO
            

C     =========== End of Main program =============

C     Now we need to put the variables in the form required for
C     output from this procedure.
      
      cbaseH2O   = cbase(1)
      cbaseNH3   = cbase(2)
      cbaseH2S   = cbase(3)  ! Not returned
      cbaseNH4SH = cbase(4)

      shH2O   = sh(1)
      shNH3   = sh(2)
      shH2S   = sh(3)   ! Not returned
      shNH4SH = sh(4)

      DO level = 1,NPRO
         vmrH2O(level) = vmr(1,level)
         vmrNH3(level) = vmr(2,level)
         vmrH2S(level) = vmr(3,level)         
      ENDDO

      END
C2345678901234567890123456789012345678901234567890123456789012345678901234567890
C        1         2         3         4         5         6         7 *
