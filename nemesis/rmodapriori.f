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
      integer varident(mvar,3),ivar,nvar,nlevel,npvar
      character*100 opfile,buffer,ipfile
 
      do i=1,nx
        err(i)=sqrt(sx(i,i))
      enddo 

      ix=0
      do 10 ivar=1,nvar
          print*,(varident(ivar,i),i=1,3)
          print*,(varparam(ivar,i),i=1,5)
          if((varident(ivar,3).eq.0).or.(varident(ivar,3).eq.25))then
C         ********* continuous profile ************************
            if(varident(ivar,3).eq.25)then
             np = varparam(ivar,1)
             clen = 1.5
            else
             np = npro
             xfac = sx(ix+1,ix+2)/sqrt(sx(ix+1,ix+1)*sx(ix+2,ix+2))
             clen = 1.0/sqrt(-alog(xfac))
            endif
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
          else
            np=1

            if(varident(ivar,1).le.100)then
             np=npvar(varident(ivar,3),npro,varparam(ivar,1))
            endif
            if(varident(ivar,1).eq.888)np=int(varparam(ivar,1))
            if(varident(ivar,1).eq.887)np=int(varparam(ivar,1))
            if(varident(ivar,1).eq.444)np=2+int(varparam(ivar,1))
            if(varident(ivar,1).eq.445)np=3+(2*int(varparam(ivar,1)))
            if(varident(ivar,1).eq.222)np = 8
            if(varident(ivar,1).eq.223)np = 9
            if(varident(ivar,1).eq.224)np = 9
            if(varident(ivar,1).eq.225)np = 11
            if(varident(ivar,1).eq.226)np = 8
            if(varident(ivar,1).eq.227)np = 7
            

            do i=1,np
             dx=2*(ran11(idum)-0.5)*err(ix+i)             
             xn(ix+i)=x0(ix+i)+dx
            enddo

          endif

          ix = ix+np

10    continue

      return

      end

