      subroutine modxvecMCMC(idum,npro,nvar,varident,varparam,
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
      real efsh,xfsh,varparam(mvar,mparam),gasdev1
      real ref(maxpro),clen,pi,dx,xl,xphi,xamp,sum,f
      parameter(pi=3.1415927)
      integer varident(mvar,3),ivar,nvar,nlevel
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
      do 10 ivar=1,nvar
          if(varident(ivar,3).eq.0)then
C         ********* continuous profile ************************

C           Randomly modify the trial profile
C            print*,ivar,npro
            np=npro
            do i=1,np
             xtry(i)=err(ix+i)*gasdev1(idum)
C             print*,i,x0(ix+i),err(ix+i),xtry(i)
            enddo

C           Now smooth with covariance matrix
            do i=1,np
             xn(ix+i)=0.
             sum=0.
             do j=1,np
               f=sqrt(sx(ix+i,ix+j)/sx(ix+i,ix+i))
C              print*,j,f
               xn(ix+i)=xn(ix+i) + f*xtry(j)
              sum=sum+f
             enddo
             xn(ix+i)=x0(ix+i)+xn(ix+i)/sum
            enddo          

      
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
          endif

          ix = ix+np

10    continue

      return

      end