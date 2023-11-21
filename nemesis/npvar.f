      integer function npvar(imod,npro,vpar1)
C     **************************************************************
C     Simple function to define the number of elements in the measurement
C     vector, xn, needed to describe an atmospheric profile, depending on
C     the require parameterisation scheme.
C
C     Input variable
C	imod	integer	Require parameterisation scheme ID.
C  	npro	integer	Number of vertical levels in .ref file
C	vpar1	real	varparam(ivar,1)
C
C     Output variable
C	npvar	integer	Number of elements in xn required.
C
C     Pat Irwin	11/5/12
C
C     **************************************************************
      implicit none
      integer imod,np,npro,nlong
      real vpar1
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      np=1
      if(imod.le.50)then
        if(imod.eq.-1)np = npro
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
        if(imod.eq.24)np = 3
        if(imod.eq.25)np = int(vpar1)
        if(imod.eq.26)np = 4
        if(imod.eq.27)np = 3
        if(imod.eq.28)np = 1
        if(imod.eq.29)np = int(vpar1)*npro
        if(imod.eq.30)np = int(vpar1)
        if(imod.eq.31)np = int(vpar1)
        if(imod.eq.32)np = 3
        if(imod.eq.33)then
          nlong=16
          if(idiag.gt.0)then
           print*,'npvar - model 33: assuming nlong = ',nlong
          endif
          np = int(vpar1)+nlong
        endif
        if(imod.eq.34)np = 2
        if(imod.eq.35)np = int(vpar1)
        if(imod.eq.36)np = int(vpar1)
        if(imod.eq.40)np = 2
        if(imod.eq.41)np = 5
        if(imod.eq.42)np = 3
        if(imod.eq.43)np = 5
        if(imod.eq.44)np = 5*int(vpar1)
        if(imod.eq.45)np = 3
        if(imod.eq.46)np = 6
        if(imod.eq.47)np = 3
        if(imod.eq.48)np = 4
        if(imod.eq.50)np = 3
      else
        print*,'npvar:  Model parameterisation not defined = ',imod
      endif

      npvar=np

      return

      end


