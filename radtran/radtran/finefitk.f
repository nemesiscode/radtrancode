      subroutine finefitk(ndata,u,trans,etrans,ng,del_g,k_g,icalc)
C     ************************************************************************
C     Fits and exponential sum series to an input transmission curve using
C     optimal estimation retrieval theory.
C
C     Input Variables:
C	ndata	     	integer	Number of points in transmission curve
C	u(mdata)	real	Amounts
C	trans(mdata)	real	Transmissions
C	etrans(mdata)	real	Transmission errors
C	ng		integer	Number of k-coefficients
C	del_g(10)	real	Weights of coefficients
C	k_g(10)		real	First guess k-coefficients
C
C     Output variables
C	k_g(10)		real	Fitted coefficients
C	icalc		integer	Quality indicator:
C				0 : successful fit
C				-1 : May not have converged
C
C     Pat Irwin
C       Original July 2003
C	Revised and documented	9/1/04
C	Revised again		24/1/04
C       Overhauled		20/8/13
C
C     ***********************************************************************
      implicit none
      include '../includes/arrdef.f'
      include '../../nemesis/arraylen.f'
      integer mdata
      parameter (mdata=20)
      real u(mdata),trans(mdata),etrans(mdata),tmp
      integer ng,ndata,kiter,i,j,iter
      real alambda,xchi,oxchi,phlimit
      real del_g(maxg),k_g(maxg),tphi
      real kk(my,mx),sa(mx,mx),se(my,my),xn(mx),xa(mx),y(my)
      double precision sei(my,my),sai(mx,mx)
      real xn1(mx),x_out(mx),yn1(my),xx
      integer nx,ny,icalc,iflag
      real phi,ophi,calc_phiret,kk1(my,mx),yn(my),chisq
      double precision aa(mx,mx),dd(mx,my)
C     -----------------------------------------------------------------------

      nx=ng
      ny=ndata
      phlimit=0. 
      icalc=0

C     Reset arrays
      do i=1,mx
       do j=1,mx
        sa(i,j)=0.
       enddo
       do j=1,my
        kk(j,i)=0.
        kk1(j,i)=0.
       enddo
      enddo

      do i=1,my
       do j=1,my
        se(i,j)=0.
        sei(i,j)=0.
       enddo
      enddo

      do 10 i=1,nx
       xa(i)=alog(k_g(i))
       xn(i)=xa(i)
C      xn is log(k), so errors are fractional errors. Set to 1.0 to give model
C      free reign in a least squares fit.
       sa(i,i)=1.
       sai(i,i)=1./sa(i,i)
10    continue

      do 11 j=1,ny
       y(j)=alog(trans(j))
       se(j,j)=1e-6
       sei(j,j)=1.0/se(j,j)
11    continue


      
C     Calculate first transmission spectrum
      call calc_trans(nx,ny,u,xn,del_g,yn,kk)

C     Now calculate the gain matrix and averaging kernels
      call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)

C     Calculate initial value of cost function phi.
      phi = calc_phiret(ny,y,yn,sei,nx,xn,xa,sai,chisq)
      ophi = phi

      print*,'Initial phi = ',phi

C     Assess whether retrieval is likely to be OK
      call assess(nx,ny,kk,sa,se)

C     alambda is a Marquardt-Levenberg-type 'braking parameter'
      alambda = 1.0

C     Set the trial vectors xn1, and yn1 to be the same as the initial
C     vectors xn, yn
      do i=1,nx
       xn1(i)=xn(i)
      enddo
      do i=1,ny
       yn1(i)=yn(i)
      enddo

      kiter=15
      do 401 iter = 1, kiter

C      Now calculate next iterated xn1
       call calcnextxn(nx,ny,xa,xn,y,yn,dd,aa,x_out)

C      Check monotonacity
901    continue
       iflag=0
       do i=1,nx-1
        if(x_out(i).gt.x_out(i+1))then
         xx=x_out(i)
         x_out(i)=x_out(i+1)
         x_out(i+1)=xx
         iflag=1
        endif
       enddo
       if(iflag.eq.1) goto 901

      do i=1,nx
         xn1(i) = xn(i) + (x_out(i)-xn(i))/(1.0+alambda)
       enddo


C      Calculate new transmission spectrum from this trial solution
       call calc_trans(nx,ny,u,xn1,del_g,yn1,kk1)
 

C      Calculate the cost function for this trial solution.
       phi = calc_phiret(ny,y,yn1,sei,nx,xn1,xa,sai,chisq)

       print*,'alambda, phi, ophi',alambda,phi,ophi

C      Does trial solution fit the data better?
       if(phi.le.ophi)then
          print*,'Successful iteration. Updating xn,yn and kk'
          do i=1,nx
           xn(i)=xn1(i)                         ! update xn to new value
          enddo
          do i=1,ny
           yn(i)=yn1(i)                         ! update yn and kk
           do j=1,nx
            kk(i,j)=kk1(i,j)
           enddo
          enddo

C         Now calculate the gain matrix and averaging kernels
          call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)      

C         Has solution converged?
          tphi = 100.0*(ophi-phi)/ophi
          if(tphi.ge.0.0.and.tphi.le.phlimit.and.alambda.lt.1.0)then
            print*,'%phi, phlimit : ',tphi,phlimit
            print*,'Phi has converged'
            print*,'Terminating retrieval'
            GOTO 202
          else
            ophi=phi
            oxchi = xchi
            alambda = alambda*0.3               ! reduce Marquardt brake
          endif
        else
C         Leave xn and kk alone and try again with more braking
          alambda = alambda*10.0                ! increase Marquardt brake
          if(alambda.gt.1e10)alambda=1e10
        endif

401   continue


202   continue
      return

      end


      subroutine calc_trans(nx,ny,u,xn,del_g,yn,kk)
C     ***************************************************************
C     Given k-coefficients, calculate transmission as function of path
C     amount and rate of change of these transmissions with each k-coefficient
C     
C     Input variables
C	nx	integer	Number of k-coefficients
C	ny	integer	Number of transmissions
C	u(ny)	real	Transmissions
C	xn(nx)	real	log of k-coefficients
C  	del_g(nx) real	Weights
C
C     Output variables
C	yn(ny)	real	Log of transmissions
C       kk(ny,nx) real	Jacobians
C
C     Pat Irwin   20/8/13
C
C     ***************************************************************
	
      implicit none
      include '../includes/arrdef.f'
      include '../../nemesis/arraylen.f'
      integer nx,ny,i,j
      real u(20),xn(mx),del_g(10),yn(my),kk(my,mx),xtrans,kx


C     Calculate transmission spectrum (as function of u)
      do 15 i=1,ny
       xtrans=0.
       do 20 j=1,nx
        kx=exp(xn(j))
        xtrans=xtrans+del_g(j)*exp(-kx*u(i))
20     continue
       yn(i)=alog(xtrans)
       do 25 j=1,nx
        kx=exp(xn(j))
        kk(i,j)= -kx*(u(i)/xtrans)*del_g(j)*exp(-kx*u(i))
25     continue
15    continue

      return
      
      end

