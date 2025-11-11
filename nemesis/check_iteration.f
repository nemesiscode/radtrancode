      subroutine check_iteration(nvar,varident,varparam,xn,npro,
     1  icheck)
C     **********************************************************
C     Subroutine to check that profiles are still valid. Called by coreretPT.
C     **********************************************************      
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer nvar,nx,varident(mvar,3),icheck,npro,i,ilog
      real xn(mx),tmp(maxpro),varparam(mvar,mparam)
      integer np,npvar,istart,logflag,imod,itype,ivar
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      icheck=0
      istart=0
      do 1000 ivar=1,nvar
       np=1
       imod=varident(ivar,3)
       itype=varident(ivar,1)
       if(itype.le.100)then
         np = npvar(imod,npro,varparam(ivar,1))
       endif
       if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
       if(varident(ivar,1).eq.887)np = int(varparam(ivar,1))
       if(varident(ivar,1).eq.444)then
        if(varparam(ivar,2).gt.0.0)then
         np = 2+int(varparam(ivar,1))
        else
         np = 3
        endif
       endif
       if(varident(ivar,1).eq.446)then
        if(varparam(ivar,2).gt.0.0)then
         np = 3+2*int(varparam(ivar,1))
        else
         np = 5
        endif
       endif
       if(varident(ivar,1).eq.445)np = 3+(2*int(varparam(ivar,1)))
       if(varident(ivar,1).eq.222)np = 8
       if(varident(ivar,1).eq.223)np = 9
       if(varident(ivar,1).eq.224)np = 9
       if(varident(ivar,1).eq.225)np = 11
       if(varident(ivar,1).eq.225)np = 8
       if(varident(ivar,1).eq.227)np = 7
       if(varident(ivar,1).eq.11)np = 2+int(varparam(ivar,1))

       do i=1,np
         tmp(i)=xn(i+istart)
       enddo

C      See if this is a temperature model. If temperatures gone negative
C      then send error flag.
       if(itype.eq.0) then 
         do i=1,np
          ilog=logflag(itype,imod,varparam(ivar,1),i)
          if(ilog.eq.0.and.tmp(i).lt.0.)icheck=1
         enddo
       endif

       if(varident(ivar,1).ge.444.and.varident(ivar,1).le.446)then
C       Check to see if the width of the size distribution has gone too
C       small, which will muck up Mie Scattering

        if(tmp(2).lt.0.01)icheck=1

       endif

       istart=istart+np

1000  continue
 
      return

      end
