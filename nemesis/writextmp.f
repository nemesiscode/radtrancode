      subroutine writextmp(runname,xlat,xlon,nvar,varident,varparam,
     1 npro,nx,xn,sx,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac)
C     $Id:
C     ************************************************************************
C     Subroutine to write out the .str file which is stripped retrieved
C     information from a .mrp file
C
C     Input variables
C       runname         character*100    Root run name.
C	xlat		real		Central latitude
C	xlon		real		Central longitude
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

      real xlat,xlon
      integer nx,npro,jlogg,jfrac
      real xn(mx),sx(mx,mx)
      integer nvar,varident(mvar,3),i,j,jsurf,jalb,jtan,jpre,jrad
      integer jxsc
      real varparam(mvar,mparam)
      character*100 runname
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      call file(runname,runname,'str')
      open(12,file=runname,status='unknown')
        write(12,*)xlat,xlon,nvar
        do i=1,nvar
         write(12,*)(varident(i,j),j=1,3)
         write(12,*)(varparam(i,j),j=1,mparam)
        enddo
        write(12,*)npro
        write(12,*)nx
        write(12,*)(xn(j),j=1,nx)
        do i=1,nx
         write(12,*)(sx(j,i),j=1,nx)
        enddo
        write(12,*)jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac
      close(12)

      return

      end
