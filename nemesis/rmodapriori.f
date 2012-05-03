      subroutine rmodapriori(idum,npro,nvar,varident,varparam,
     1  nx,x0,sx,xn)
C     $Id:
C     ****************************************************************
C     Subroutine to randomly modify x-vector for program Generatespx
C
C     Input variables
C	idum		integer		Seed for random variables
C	opfile		character*100	root name of run.
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
C	nx 		integer 	number of elements in state vector
C	x0(mx)		real		a priori vector
C	sx(mx,mx)	real		Covariance matrix
C
C     Output variable
C	xn(mx)		real		Modified vector
C
C 
C     Original:	Pat Irwin		3/3/04
C
C     ****************************************************************

      implicit none

      integer i,j,nx,ix,jx,npro,idum,np

C     ****************************************************************
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
C     ****************************************************************

      real x0(mx),xn(mx),sx(mx,mx),pref(maxpro),eref(maxpro)
      real delp,xfac,pknee,edeep,xdeep,err(mx)
      real efsh,xfsh,varparam(mvar,mparam),ran11
      real ref(maxpro),clen,pi,dx,xl,xphi,xamp
      parameter(pi=3.1415927)
      integer varident(mvar,3),ivar,nvar,nlevel
      character*100 opfile,buffer,ipfile
 
      do i=1,nx
        err(i)=sqrt(sx(i,i))
      enddo 

      ix=0
      do 10 ivar=1,nvar
          if(varident(ivar,3).eq.0)then
C         ********* continuous profile ************************
            np = npro
            xfac = sx(ix+1,ix+2)/sqrt(sx(ix+1,ix+1)*sx(ix+2,ix+2))
            clen = 1.0/sqrt(-alog(xfac))
            print*,'Assumed correlation = ',clen
            xl = clen*(0.5 + 19.5*ran11(idum))
            xphi = 2*pi*ran11(idum)
            print*,'xl,xphi = ',xl,xphi
            xamp = ran11(idum)
            print*,'xl,xphi,xamp = ',xl,xphi,xamp
C            xamp = 2.0*(ran11(idum)-0.5)
            do i=1,np
             dx=xamp*sin(2*pi*i/xl + xphi)*err(ix+i)
C             dx=xamp*err(ix+i)
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.1)then
C         ******** profile held as deep amount, fsh and knee pressure **
            np = 2
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.2)then  
C         **** Simple scaling factor of reference profile *******
            np = 1
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.3)then  
C         **** Simple log scaling factor of reference profile *******
            np = 1
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.4)then
C         ******** profile held as deep amount, fsh and variable knee press
            np = 3
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.6)then
C         ******** profile cloud od at specified level plus fsh.
            np = 2
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.666)then  
C         **** Pressure at given altitude *******
            np = 1
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.777)then  
C         **** Tangent height correction *******
            np = 1
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.999)then  
C         **** Surface temperature *******
            np = 1
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          elseif(varident(ivar,3).eq.555)then  
C         **** Radius of planet *******
            np = 1
            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo
          endif

          ix = ix+np

10    continue

      return

      end
