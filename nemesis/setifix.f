      subroutine setifix(xa,sa,nvar,varident,varparam,npro,ifix)
C     ***************************************************************
C     Subroutine to see if the fractional error on any of the state vector
C     variables is so small that it isn't worth bothering to calculate the 
C     associated row in the kk matrix. 
C
C     Not actually applicable to non-scattering calculations (which
C     use implicit differentiation) but can usefully speed up scattering
C     retrievals if we want to fix any of the variables.
C
C     Input variables
C	xa(mx)			real	a priori state vector
C	sa(mx,mx)		real	a priori covariance matrix
C	varident(mvar,3)	integer	variable parameterisation ID
C	varparam(mvar,mparam)	real	associated variables
C	npro			integer	Number of levels in atmosphere
C
C     Output variables
C	ifix(mx)		integer	elements set to 1 if no row of the
C					kk matrix is required. Set to 0
C					otherwise.
C
C     Pat Irwin    3/7/14 
C
C     ***************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer nvar,varident(mvar,3),ifix(mx),nxtemp,ivar,i,np
      integer npvar,iflag,logflag,npro,ix,j
      real xa1,ea1,ferr,minferr
C     Set minimum fractional error to fix variable.
      parameter (minferr = 1e-6)
      real varparam(mvar,mparam),sa(mx,mx),xa(mx)

            
      nxtemp=0
      print*,'Calling setifix'
      do 299 ivar=1,nvar
       np=1
       if(varident(ivar,1).le.100)then
         np = npvar(varident(ivar,3),npro)
       endif
       if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
       if(varident(ivar,1).eq.887)np = int(varparam(ivar,1))
       if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))
       if(varident(ivar,1).eq.222)np = 8
       if(varident(ivar,1).eq.223)np = 9
       if(varident(ivar,1).eq.224)np = 9
       if(varident(ivar,1).eq.225)np = 11
       if(varident(ivar,1).eq.226)np = 8

       do i=1,np
        ix = nxtemp+i
        xa1 = xa(ix)
        ea1 = sqrt(abs(sa(ix,ix)))

        iflag = logflag(varident(ivar,1),varident(ivar,3),i)

        if(iflag.eq.1)then
          xa1 = exp(xa1)
          ea1 = xa1*ea1
        endif
   
        ferr = abs(ea1/xa1)
 
        if(ferr.le.minferr)then
         ifix(ix) = 1
        else
         ifix(ix) = 0
        endif
       enddo
       
       nxtemp=nxtemp+np

299   continue

      print*,'ix,ifix : ',ix,(ifix(j),j=1,nxtemp)

      return

      end
