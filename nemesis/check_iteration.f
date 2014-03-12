      subroutine check_iteration(nvar,varident,varparam,xn,npro,icheck)
C     **********************************************************
C     Subroutine to check that profiles are still valid. Called by coreretPT.
C     **********************************************************      
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer nvar,nx,varident(mvar,3),icheck
      real xn(mx),tmp(maxpro),varparam(mvar,mparam)
      integer np,npvar,istart,logflag

      icheck=0
      istart=0
      do 1000 ivar=1,nvar
       np=1
       imod=varident(ivar,3)
       itype=varident(ivar,1)
       if(itype.le.100)then
         np = npvar(imod,npro)
       endif
       if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
       if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))

       do i=1,np
         tmp(i)=xn(i+istart)
       enddo

C      See if this is a temperature model. If temperatures gone negative
C      then send error flag.
       if(itype.eq.0) then 
         do i=1,np
          ilog=logflag(itype,imod,ip)
          if(ilog.eq.0.and.tmp(i).lt.0.)icheck=1
         enddo
       endif
       istart=istart+np
1000  continue
 
      return

      end
