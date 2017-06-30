      integer function npvar(imod,npro)
C     **************************************************************
C     Simple function to define the number of elements in the measurement
C     vector, xn, needed to describe an atmospheric profile, depending on
C     the require parameterisation scheme.
C
C     Input variable
C	imod	integer	Require parameterisation scheme ID.
C  	npro	integer	Number of vertical levels in .ref file
C
C     Output variable
C	npvar	integer	Number of elements in xn required.
C
C     Pat Irwin	11/5/12
C
C     **************************************************************
      implicit none
      integer imod,np,npro

      np=1
      if(imod.le.23)then
        if(imod.eq.0)np = npro
        if(imod.eq.1)np = 2
        if(imod.eq.4)np = 3
        if(imod.eq.6)np = 2
        if(imod.eq.7)np = 2
        if(imod.eq.8)np = 3
        if(imod.eq.9)np = 3
        if(imod.eq.10)np = 4
        if(imod.eq.11)np = 2
        if(imod.eq.12)np = 3
        if(imod.eq.13)np = 3
        if(imod.eq.14)np = 3
        if(imod.eq.15)np = 3
        if(imod.eq.16)np = 4
        if(imod.eq.17)np = 2
        if(imod.eq.18)np = 2
        if(imod.eq.19)np = 4
        if(imod.eq.20)np = 2
        if(imod.eq.21)np = 2
        if(imod.eq.22)np = 5
        if(imod.eq.23)np = 4
      else
       print*,'npvar:  Model parameterisation not defined = ',imod
      endif

      npvar=np

      return

      end


