      subroutine readxtmp(runname,xlat,nvar,varident,varparam,npro,
     1 nx,xn,sx,jsurf,jalb,jtan,jpre,jrad,jlogg)
C     $Id:
C     ************************************************************************
C     Subroutine to read in a .str file which is stripped retrieved
C     information from a .mrp file
C
C     Input variables
C       runname         character*100    Root run name.
C     Output variables
C	xlat		real		Central latitude
C	nvar		integer		Number of variables
C       varident(mvar,3) integer 	identity of constituent to retrieved
C					 and parameterisation
C       varparam(mvar,mparam) real 	Additional parameters constraining
C					  profile.
C	nx		integer		Number of elements in state vector
C	xn(mx)		real		state vector
C	sx(mx,mx)	real		state covariance matrix
C
C     Pat Irwin	19/8/04	New
C
C     ************************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'

      real xlat
      integer nx,jlogg
      real xn(mx),sx(mx,mx)
      integer nvar,npro,varident(mvar,3),i,j,jsurf,jalb,jtan,jpre,jrad
      real varparam(mvar,mparam)
      character*100 runname

      call file(runname,runname,'str')
      open(12,file=runname,status='old')
        read(12,*)xlat,nvar
        do i=1,nvar
         read(12,*)(varident(i,j),j=1,3)
         read(12,*)(varparam(i,j),j=1,mparam)
        enddo
        read(12,*)npro
        read(12,*)nx
        read(12,*)(xn(j),j=1,nx)
        do i=1,nx
         read(12,*)(sx(j,i),j=1,nx)
        enddo
        read(12,*)jsurf,jalb,jtan,jpre,jrad,jlogg
      close(12)

      return

      end
