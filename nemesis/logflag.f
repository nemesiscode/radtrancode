      integer function logflag(ivar,imod,vpar1,ip)
C     **************************************************************
C     Simple function to define whether a particular variable is held as
C     a log variable within Nemesis.
C
C     Input variable
C	ivar	integer	Profile ID (varident(1)).
C	imod	integer	Required parameterisation scheme ID.
C       vpar1	real	First element of varparam
C	ip	integer	Element of parameterisation vector
C
C     Output variable
C	logflag	integer	Set to 1 if log, 0 otherwise.
C
C     Pat Irwin	11/5/12
C
C     **************************************************************
      implicit none
      integer ivar,imod,iflag,ip
      real vpar1
      iflag=0

      if(ivar.ne.0)then
C       Variable is not temperature  - may need to take exponent
        if(imod.eq.0)iflag=1		 ! continuous profile
        if(imod.eq.1.and.ip.eq.1)iflag=1 ! knee profile         
        if(imod.eq.20.and.ip.eq.1)iflag=1 ! knee profile         
        if(imod.eq.4.and.ip.eq.1)iflag=1 ! variable knee profile
        if(imod.eq.6.and.ip.eq.1)iflag=1 ! Venus cloud profile
        if(imod.eq.7.and.ip.eq.1)iflag=1 ! extended profile         
        if(imod.eq.17.and.ip.eq.1)iflag=1 ! extended profile         
        if(imod.eq.18.and.ip.eq.1)iflag=1 ! extended profile         
        if(imod.eq.8.and.ip.eq.1)iflag=1 ! variable knee profile
        if(imod.eq.9.and.ip.eq.1)iflag=1 ! variable knee profile
        if(imod.eq.32.and.ip.eq.1)iflag=1 ! variable knee profile
        if(imod.eq.21.and.ip.eq.1)iflag=1 ! variable knee profile
        if(imod.eq.24.and.ip.eq.1)iflag=1 ! deep profile
        if(imod.eq.27.and.ip.eq.1)iflag=1 ! step profile (deep value)
        if(imod.eq.27.and.ip.eq.2)iflag=1 ! step profile (shallow value)
        if(imod.eq.16.and.ip.eq.1)iflag=1 ! lapse rate profile      
        if(imod.eq.19.and.ip.eq.1)iflag=1 ! lapse rate profile      
        if(imod.eq.25)iflag=1  ! Shortened continuous model
        if(imod.eq.28)iflag=1 !Modify just one element of a profile
        if(imod.eq.29)iflag=1 ! continuous profile
        if(imod.eq.30)iflag=1 ! continuous profile inhomogeneous disc
      endif

      if(imod.eq.31)iflag=1 ! log multiplier - inhomogeneous disc
      if(imod.eq.1.and.ip.eq.2)iflag=1 ! log fsh - fixed knee
      if(imod.eq.20.and.ip.eq.2)iflag=1 ! log fsh - fixed knee
      if(imod.eq.4.and.ip.eq.2)iflag=1 ! log fsh - var. knee
      if(imod.eq.4.and.ip.eq.3)iflag=1 ! variable knee profile
      if(imod.eq.6.and.ip.eq.2)iflag=1 ! Venus cloud profile
      if(imod.eq.7.and.ip.eq.2)iflag=1 ! log fsh - extended
      if(imod.eq.17.and.ip.eq.2)iflag=1 ! log fsh - extended
      if(imod.eq.18.and.ip.eq.2)iflag=1 ! log fsh - extended
      if(imod.eq.8.and.ip.eq.2)iflag=1 ! log fsh - var. knee
      if(imod.eq.8.and.ip.eq.3)iflag=1 ! variable knee profile
      if(imod.eq.32.and.ip.eq.2)iflag=1 ! variable knee profile
      if(imod.eq.32.and.ip.eq.3)iflag=1 ! variable knee profile
      if(imod.eq.9.and.ip.eq.2)iflag=1 ! log fsh - var. knee
      if(imod.eq.19.and.ip.eq.2)iflag=1 ! log fsh - var. knee
      if(imod.eq.9.and.ip.eq.4)iflag=1 ! log cwid - var knee
      if(imod.ge.12.and.imod.le.13)then 
          if(ip.eq.1)iflag=1 	       ! Gaussian/Lorentz cloud
          if(ip.eq.2)iflag=1 
          if(ip.eq.3)iflag=1 
      endif
      if(imod.ge.14.and.imod.le.15)then 
          if(ip.eq.1)iflag=1 	       ! Gaussian/Lorentz cloud
          if(ip.eq.3)iflag=1 
      endif
      if(imod.eq.16)then
          if(ip.eq.2)iflag=1 	       ! Lapse rate profile
          if(ip.eq.3)iflag=1 
          if(ip.eq.4)iflag=1 
      endif

      if(imod.eq.24.and.ip.eq.2)iflag=1 ! profile between knee and condensation
      if(imod.eq.24.and.ip.eq.3)iflag=1 ! knee pressure
      
      if(imod.eq.27.and.ip.eq.3)iflag=1 ! knee pressure for step profile
      
      if(imod.eq.-1)iflag=1 ! Dust continuous profile in particles/cm3

      if(imod.eq.3)iflag=1	! Log scaling factor
      if(imod.eq.10)iflag=1	! Log scaling factor
      if(imod.eq.11)iflag=1	! Log scaling factor

      if(imod.eq.22)then 
       iflag=1			! Brown dwarf T-profile
      endif
      
      if(imod.eq.23)iflag=1  ! 2 point vmr gradient
      if(imod.eq.26)iflag=1  ! 2 point vmr gradient (zero vmr for deep atmos)

      if(ivar.eq.887)iflag=1	! X-section spectrum
      if(ivar.eq.888)iflag=1	! Surface albedo spectrum
      if(ivar.eq.889)iflag=1	! Surface albedo spectrum multiplier
      if(ivar.eq.666)iflag=1	! Tangent pressure
      if(ivar.eq.555)iflag=0	! Planet radius
      if(ivar.eq.333)iflag=0	! Planet surface gravity
      if(ivar.eq.444)iflag=1	! Particle size and ref. index
      if(ivar.eq.445)iflag=1	! Particle size and ref. index (coated sphere)
      if(ivar.eq.443)then
         if(ip.eq.1)iflag=1     ! Log scaling factor
         if(ip.eq.2)iflag=0     ! Knee height
         if(ip.eq.3)iflag=0     ! Power law exponent
      endif
      if(ivar.eq.442)then
         if(ip.eq.1)iflag=1     ! Log scaling factor
         if(ip.eq.2)iflag=0     ! Knee height
         if(ip.eq.3)iflag=0     ! Knee height
         if(ip.eq.4)iflag=0     ! Power law exponent
      endif
      if(ivar.eq.441)then
         if(ip.eq.1)iflag=1     ! Log scaling factor
         if(ip.eq.2)iflag=0     ! Knee height
         if(ip.eq.3)iflag=0     ! Power law exponent
      endif
      if(ivar.eq.222)iflag=1	! Larry's cloud model
      if(ivar.eq.223)iflag=1	! Larry's revised cloud model
      if(ivar.eq.224)iflag=1	! Larry's revised cloud model with ext UTC
      if(ivar.eq.225)iflag=1	! Revised cloud model with ext UTC and trunk.
      if(ivar.eq.226)iflag=1	! Two cloud model
      if(ivar.eq.227)iflag=1	! Creme Brulee
      logflag=iflag

      return

      end


