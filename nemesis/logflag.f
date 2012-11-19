      integer function logflag(ivar,imod,ip)
C     **************************************************************
C     Simple function to define whether a particular variable is held as
C     a log variable within Nemesis.
C
C     Input variable
C	ivar	integer	Profile ID (varident(1)).
C	imod	integer	Required parameterisation scheme ID.
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

      iflag=0

      if(ivar.ne.0)then
C       Variable is not temperature  - may need to take exponent
        if(imod.eq.0)iflag=1		 ! continuous profile
        if(imod.eq.1.and.ip.eq.1)iflag=1 ! knee profile         
        if(imod.eq.4.and.ip.eq.1)iflag=1 ! variable knee profile
        if(imod.eq.6.and.ip.eq.1)iflag=1 ! Venus cloud profile
        if(imod.eq.7.and.ip.eq.1)iflag=1 ! extended profile         
        if(imod.eq.8.and.ip.eq.1)iflag=1 ! variable knee profile
        if(imod.eq.9.and.ip.eq.1)iflag=1 ! variable knee profile
      endif

      if(imod.eq.1.and.ip.eq.2)iflag=1 ! log fsh - fixed knee
      if(imod.eq.4.and.ip.eq.2)iflag=1 ! log fsh - var. knee
      if(imod.eq.4.and.ip.eq.3)iflag=1 ! variable knee profile
      if(imod.eq.6.and.ip.eq.2)iflag=1 ! Venus cloud profile
      if(imod.eq.7.and.ip.eq.2)iflag=1 ! log fsh - extended
      if(imod.eq.8.and.ip.eq.2)iflag=1 ! log fsh - var. knee
      if(imod.eq.8.and.ip.eq.3)iflag=1 ! variable knee profile
      if(imod.eq.9.and.ip.eq.2)iflag=1 ! log fsh - var. knee
      if(imod.ge.12.and.imod.le.13)then 
          if(ip.eq.1)iflag=1 	       ! Gaussian/Lorentz cloud
          if(ip.eq.2)iflag=1 
          if(ip.eq.3)iflag=1 
      endif
      if(imod.ge.14.and.imod.le.15)then 
          if(ip.eq.1)iflag=1 	       ! Gaussian/Lorentz cloud
          if(ip.eq.3)iflag=1 
      endif

      if(imod.eq.3)iflag=1	! Log scaling factor
      if(imod.eq.10)iflag=1	! Log scaling factor
      if(imod.eq.11)iflag=1	! Log scaling factor

      if(ivar.eq.888)iflag=1	! Surface albedo spectrum
      if(ivar.eq.666)iflag=1	! Tangent pressure

      logflag=iflag

      return

      end


