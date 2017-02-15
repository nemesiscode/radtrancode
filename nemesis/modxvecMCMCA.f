      subroutine modxvecMCMCA(idum,npro,nvar,varident,varparam,
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
C	x0(mx)		real		Previous measurement vector
C	sx(mx,mx)	real		Covariance of proposal distribution
C
C     Output variable
C	xn(mx)		real		Modified vector
C
C 
C     Original:	Pat Irwin		3/3/04
C     Modifed from rmodaprioiri.f Pat Irwin 11/6/13
C
C     ****************************************************************

      implicit none

      integer i,j,nx,ix,jx,npro,idum,np

C     ****************************************************************
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
C     ****************************************************************

      real x0(mx),xn(mx),sx(mx,mx),pref(maxpro),eref(maxpro)
      real delp,xfac,pknee,edeep,xdeep,err(mx),xtry(mx)
      real efsh,xfsh,varparam(mvar,mparam),gasdev1,emax
      real ref(maxpro),clen,pi,dx,xl,xphi,xamp,sum,f
      real xtry2(mx),xx(mx),xy,xnow
      parameter(pi=3.1415927)
      integer varident(mvar,3),ivar,nvar,nlevel,nspline
      character*100 opfile,buffer,ipfile


C      print*,idum,npro,nvar
C      do i=1,nvar
C       print*,(varident(i,j),j=1,3)
C       print*,(varparam(i,j),j=1,5)
C      enddo
C      print*,nx
C      print*,(x0(i),i=1,nx)

      do i=1,nx
        err(i)=sqrt(sx(i,i))
      enddo
 
      ix=0
      nspline=10
      do 10 ivar=1,nvar
          if(varident(ivar,3).eq.0)then
C         ********* continuous profile ************************

C           Randomly modify the trial profile
C            print*,ivar,npro
            np=npro
            emax=0.
            do i=1,nspline
             xx(i)=100.*float(i-1)/float(nspline-1)
             j = int(1 + (np-1)*float(i-1)/float(nspline-1))
             xtry(i)=err(ix+j)*gasdev1(idum)
            enddo
            call cspline(xx,xtry,nspline,1.e30,1.e30,xtry2)

            do i=1,np
             xnow = 100*float(i-1)/float(np-1)
             call csplint(xx,xtry,xtry2,nspline,xnow,xy)
             xn(ix+i)=x0(ix+i)+xy  
            enddo          

            write(41,*)(xn(ix+j)-x0(ix+j),j=1,np)

      
          elseif(varident(ivar,3).eq.1)then
C         ******** profile held as deep amount, fsh and knee pressure **
            np = 2
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.2)then  
C         **** Simple scaling factor of reference profile *******
            np = 1
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.3)then  
C         **** Simple log scaling factor of reference profile *******
            np = 1
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.4)then
C         ******** profile held as deep amount, fsh and variable knee press
            np = 3
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.6)then
C         ******** profile cloud od at specified level plus fsh.
            np = 2
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.666)then  
C         **** Pressure at given altitude *******
            np = 1
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.777)then  
C         **** Tangent height correction *******
            np = 1
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.999)then  
C         **** Surface temperature *******
            np = 1
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.555)then  
C         **** Radius of planet *******
            np = 1
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.333)then  
C         **** surface gravity of planet *******
            np = 1
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.222)then  
C         **** Larry's cloud model *******
            np = 8
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.223)then  
C         **** Larry's revised cloud model *******
            np = 9
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.224)then  
C         **** Larry's revised cloud model with extended UTC *******
            np = 9
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.225)then  
C         **** Larry's revised cloud model with extended UTC and trunk *******
            np = 11
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.226)then  
C         **** Two cloud model  *******
            np = 8
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.888)then
C         **** Surface albedo spectrum of planet *******
            np = int(varparam(ivar,1))
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          elseif(varident(ivar,3).eq.444)then
C         **** Particle refractive index spectrum *******
            np = 2+int(varparam(ivar,1))
            do i=1,np
             xn(ix+i)=x0(ix+i)+err(ix+i)*gasdev1(idum)
            enddo
          endif

          ix = ix+np

10    continue

      return

      end
