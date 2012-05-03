      subroutine map2xvec(nparam,npro,npath,nv,xmap,doutmoddq,doutdx)
C     ****************************************************************
C     Subroutine to convert from rate of change of radiance with profile
C     .prf properties to user defined variables.
C
C     Input variables
C 	nparam	integer	Number of profile gradients defined.
C	npro	integer	Number of profile levels
C	npath	integer	Number of calculated paths
C	nv	integer	Number of user-defined variables
C	xmap	real	Matrix to convert form profile gradients to 
C			user defined gradients. 
C       doutmoddq(maxpat,maxgas+2+maxcon,maxpro) real Rate of change
C                                       of output with respect to .prf
C                                       profile properties:
C
C       1 to ngas                       profile gas vmrs
C       ngas+1                          profile temperatures
C       ngas+2 to ngas+2+ncont-1        profile cloud specific densities
C       ngas+ncont+2                    profile para-H2 fraction.
C
C     Output variable
C 	doutdx(maxpat,maxv) real	Rate of change of output with 
C					respect to each user-defined
C					variable.
C
C     Pat Irwin	30/7/01	Original version
C     Pat Irwin	29/2/12	Updated for Radtrans2.0

C     ****************************************************************

      implicit none
      include '../includes/arrdef.f'
    
      integer nparam,npro,npath,nv
      real xmap(maxv,maxgas+2+maxcon,maxpro) 
      real doutmoddq(maxpat,maxgas+2+maxcon,maxpro)
      real doutdx(maxpat,maxv)
      integer ipath,iv,iparam,ipro

      do 100 ipath=1,npath
        do 10 iv=1,nv

	 doutdx(ipath,iv)=0.0

	 do iparam=1,nparam
          do ipro=1,npro
           doutdx(ipath,iv)=doutdx(ipath,iv) +
     1	     doutmoddq(ipath,iparam,ipro)*xmap(iv,iparam,ipro)
          enddo
         enddo

10      continue
100   continue

      return

      end
