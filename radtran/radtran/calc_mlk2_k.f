      subroutine calc_mlk2_k(SL,BL,u,g_ord,k_g,ng)
C     ************************************************************************
C
C     Calculates the mean value of k(g) in a g-interval using the method
C     described in Lacis and Oinas (LO), p9038. 
C
C     g(k) and h(k) are determined for each ordinate and then kbar 
C     calculated from:
C		    1
C		k = - log(delta_g/delta_h)
C		    u
C
C	Pat Irwin	26/10/94
C
C     **********************************************************************
      implicit none
      integer maxg,i,niter,j,c1,c2,maxint,ng
      parameter (maxg=21,maxint=200)
      real g_ord(maxg),delg(maxg),u,k_g(maxg)
      real dh,tran_ml,calc_mlh,calc_mlg,tranu,k,calc_mlk0
      real lkr(maxint),lgr(maxint),kbar,kbar0
      double precision aa,bb,SL,BL
      real kmax,dg,gg,kn,kn1,sum,dlk,frac,gk,grad,pi,k1,k2,k3
      real f1,f2,f3,slvf,ke1,ke2,lk1,lk2,lk
      real kr(maxg),gr(maxg),ka(maxg),hr(maxg)
      parameter (pi=3.1415927)
      logical flag

      do i=1,ng
       gr(i)=g_ord(i)
      end do
 
      do i=1,ng-1
       delg(i)=gr(i+1)-gr(i)
      end do

C     Determine 'a' and 'b' as defined on p9037 of LO

      IF(SL.EQ.0.OR.BL.EQ.0.)THEN
       do i=1,ng-1
        k_g(i)=0.
       end do
       return
      END IF


      aa = 0.5*sqrt(pi*BL*SL)
      bb = 0.5*sqrt(pi*BL/SL)


C     Calculate the maximum value of k in the frequency distribution

      kmax = sngl((3.*SL/(pi*BL))*(sqrt((pi*BL/3.)**2 + 1.) - 1.))
     

      if(kmax.eq.0.)then
       print*,'Setting k(g) to zero'
       do i=1,ng-1
        k_g(i)=0.
       end do
       return
      end if


C     First, roughly divide up the interval
      

      ke1=kmax
      ke2=kmax
      c1=0
      c2=0
      
      gg=0.001

C     The function slvf calculates the function (g(ke1) - gg) given the LO
C      parameters 'a', 'b' and 'B' 

13    if (slvf(ke1,aa,bb,BL,gg).gt.0) then
       ke1=ke1*0.5
       c1=c1+1
       if(c1.lt.500)then
        goto 13
       else
        print*,'Error. ke1 has been halved more than 500 times'
        print*,'ke1,aa,bb,BL,gg'
        print*,ke1,aa,bb,BL,gg
        stop
       end if
      end if


      gg=0.999

23    if (slvf(ke2,aa,bb,BL,gg).lt.0) then
        ke2=ke2*2.
        c2=c2+1 
        if(c2.lt.500)then
         goto 23
        else
         print*,'Error. ke2 has been doubled more than 500 times'
         print*,'ke2,aa,bb,BL,gg'
         print*,ke2,aa,bb,BL,gg
         stop
       end if
      end if
     
      gg=0.

C     split up range ke1 to ke2 into maxint intervals and determine (maxint
C     +1) ordinates and determine for each ordinate k, g(k)
C     Values go into lkr(i),lgr(i) and area(i)

      dlk=(ke2/ke1)**(1./real(maxint-1))

      lkr(1)=ke1
      
      do 54 i=2,maxint
       k1=lkr(i-1)*dlk
       lkr(i)=k1
       lgr(i)=slvf(k1,aa,bb,BL,gg)
54    continue

C     forcing extreme limits of g(k) to be 0. and 1.
      lgr(1)=0.
      lgr(maxint)=1.
         

C     calculate k(g)

      kr(1)=ke1
      kr(ng)=ke2

      do 11 i=2,ng-1
         do 77 j=1,maxint
          if(lgr(j).ge.gr(i))then
           frac=(gr(i)-lgr(j-1))/(lgr(j) - lgr(j-1))
           kr(i)=lkr(j-1)+frac*(lkr(j)-lkr(j-1))
           goto 78
          end if
77       continue
78     continue
11    continue


      if(SL.eq.0.or.BL.eq.0.or.u.eq.0.)then
       tranu=1.
      else
       tranu=tran_ml(SL,BL,u)
      end if


      if(tranu.eq.1.)then
       print*,'TRAN_ML(SL,BL,u)=1. Hence meank = 0.'
       do i=1,ng-1
        k_g(i)=0.
       end do
       return
      endif

      if(tranu.eq.0.)then
       print*,'TRAN(SL,BL,u)=0. Hence set meank to be average'
       do i=1,ng-1
        k_g(i)=0.5*(kr(i)+kr(i+1))
       end do
       return
      end if

      do 12 i=1,ng
       k=kr(i)
       gr(i)=calc_mlg(k,SL,BL)
       hr(i)=calc_mlh(k,SL,BL,u)
12    continue

      do 14 i=1,ng-1
       dg=gr(i+1)-gr(i)
       dh=hr(i+1)-hr(i)
       k1=kr(i)
       k2=kr(i+1)
       kbar0=calc_mlk0(k1,k2,SL,BL)
       if(dh.gt.0.)then
        kbar=(1/u)*log(dg/dh)
       else
        kbar=kbar0
       end if
       if(kbar.lt.kr(i))kbar=kbar0
       if(kbar.gt.kr(i+1))kbar=kbar0
       k_g(i)=kbar
14    continue

      return

      end

      real function calc_mlg(k,SL,BL)
C     $Id: calc_mlk2_k.f,v 1.2 2011-06-17 15:40:25 irwin Exp $
C     ***********************************************************************
C
C     Calculates g(k) for a Malkmus-Lorentz band according to Lacis and Oinas 
C     conventions.
C
C     Pat Irwin		26/10/94
C
C     ***********************************************************************
      implicit none
      real k,pi,erf
      parameter (pi=3.1415927)
      double precision aa,bb,SL,BL
      real xx,yy

      aa=0.5*dsqrt(pi*BL*SL)
      bb=0.5*dsqrt(pi*BL/SL)

      xx=sngl(aa/sqrt(k))
      yy=sngl(bb*sqrt(k))

      calc_mlg=0.5*(1.-erf(xx-yy)) +0.5*(1-erf(xx+yy))*exp(pi*sngl(BL))

      return

      end


      real function calc_mlh(k,SL,BL,u)
C     $Id: calc_mlk2_k.f,v 1.2 2011-06-17 15:40:25 irwin Exp $
C     *******************************************************************
C
C     Calculates h(k) for a Malkmus-Lorentz band according to Lacis and Oinas 
C     conventions.
C
C     Pat Irwin		26/10/94
C
C     *******************************************************************
      implicit none
      real k,pi,erf,tranu,tran_ml,u
      parameter (pi=3.1415927)
      double precision aa,cc,SL,BL
      real xx,yy


      aa=0.5*dsqrt(pi*BL*SL)
      cc=0.5*dsqrt((pi*BL/SL) + 4.*u)

      xx=sngl(aa/sqrt(k))
      yy=sngl(cc*sqrt(k))

      tranu=tran_ml(SL,BL,u)

      calc_mlh=0.5*(1.-erf(xx-yy))*tranu +0.5*(1-erf(xx+yy))*
     1exp(pi*sngl(BL))/tranu

      return

      end

      real function calc_mlk0(k1,k2,SL,BL)
C     $Id: calc_mlk2_k.f,v 1.2 2011-06-17 15:40:25 irwin Exp $
C     ***********************************************************************
C
C     Calculates the mean value of k over the interval k1 to k2 using the
C     Malkmus-Lorentz band (Lacis and Oinas form) in the limit where the
C     absorber amount tends to zero. (LO p9038)
C
C     Pat Irwin		26/10/94
C
C     ***********************************************************************
      implicit none
      real k1,k2,pi,erf,tranu
      parameter (pi=3.1415927)
      double precision aa,bb,SL,BL
      real xx,yy,A1,B1


      aa=0.5*dsqrt(pi*BL*SL)
      bb=0.5*dsqrt(pi*BL/SL)

      xx=sngl(aa/sqrt(k2)) - sngl(bb*sqrt(k2))
      yy=sngl(aa/sqrt(k2)) + sngl(bb*sqrt(k2))

      A1=erf(xx)
      B1=erf(yy)


      xx=sngl(aa/sqrt(k1)) - sngl(bb*sqrt(k1))
      yy=sngl(aa/sqrt(k1)) + sngl(bb*sqrt(k1))

      A1=A1-erf(xx)
      B1=B1-erf(yy)


      calc_mlk0=sngl(SL)*(A1 - B1*exp(pi*sngl(BL)))/
     1(A1 + B1*exp(pi*sngl(BL)))

      return

      end


