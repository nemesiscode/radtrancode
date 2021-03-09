      subroutine writemvec(ivec,npro,nvar,varident,varparam,jsurf,
     1  nx,x0,xn)
C     $Id:
C     ****************************************************************
C     Subroutine to write out x-vector and associated parameters for 
C     program Randomspec
C
C     Input variables
C	ivec		integer		File ID
C	npro		integer		Number if vertical levels in .prf
C	nvar		integer		Number of variable profiles
C					(including T, vmr and cloud)
C	varident(mvar,3)integer		identity of constituent to 
C  					retrieved and how it is represented
C					First and second column contains
C					identity. Third column contains:
C					0 read in new profile and error
C					1 read in deep, fsh, knee
C					2 scale profile in .ref file.
C       varparam(mvar,mparam) integer   Additional parameters constraining
C					 profile.
C	jsurf		integer		Position of surface temperature 
C					element (if included)
C	nx 		integer 	number of elements in state vector
C	x0(mx)		real		a priori state vector
C	xn(mx)		real		actual state vector
C
C 
C     Original:	Pat Irwin		3/3/04
C
C     ****************************************************************

      implicit none

      integer i,j,nx,ix,jx,npro,jsurf,ivec

C     ****************************************************************
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
C     ****************************************************************

      real xn(mx),err,pref(maxpro),eref(maxpro)
      real delp,xfac,pknee,edeep,xdeep,x0(mx)
      real efsh,xfsh,varparam(mvar,mparam)
      real ref(maxpro),clen
      integer varident(mvar,3),ivar,nvar,nlevel
      character*100 opfile,buffer,ipfile
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet
 
      write(ivec,*)nx,'    ! nx'
      write(ivec,*)jsurf,'    ! jsurf'

      do 20 i=1,nx
       write(ivec,*)x0(i),xn(i)
20    continue

      return

      end
