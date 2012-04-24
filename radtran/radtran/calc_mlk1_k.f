      subroutine calc_mlk1_k(SL,BL,g_ord,k_g,ng)
C     ************************************************************************
C
C
C     ************************************************************************
      implicit none
      include '../includes/arrdef.f'
      integer i,niter,j,c1,c2,maxint,igdist,ng
      parameter (maxint=200)
      real g_ord(maxg),delg(maxg),area(maxint),k_g(maxg)
      real lkr(maxint),lgr(maxint)
      double precision aa,bb,SL,BL
      real kmax,dg,gg,kn,kn1,sum,dlk,frac,gk,grad,pi,k1,k2,k3
      real f1,f2,f3,slvf,ke1,ke2,lk1,lk2,lk
      real kr(maxg),gr(maxg),ka(maxg)
      parameter (pi=3.1415927)
      logical flag


      do i=1,ng
       gr(i)=g_ord(i)
      end do

      do i=1,ng-1
       delg(i)=gr(i+1)-gr(i)
      end do   

      IF(SL.EQ.0.OR.BL.EQ.0.)THEN
       do i=1,ng-1
        k_g(i)=0.
       end do
       return
      END IF

C     Determine 'a' and 'b' as defined on p9037 of LO

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
C     +1) ordinates and determine for each ordinate k, g(k) and the cumulative
C     area under g(k). Values go into lkr(i),lgr(i) and area(i)

      dlk=(ke2/ke1)**(1./real(maxint-1))

      lkr(1)=ke1
      area(1)=0.
      
      do 54 i=2,maxint
       k1=lkr(i-1)*dlk
       lkr(i)=k1
       lgr(i)=slvf(k1,aa,bb,BL,gg)
       area(i)=area(i-1) + 0.5*(lkr(i)+lkr(i-1))*(lgr(i)-lgr(i-1))
54    continue

C     forcing extreme limits of g(k) to be 0. and 1.
      lgr(1)=0.
      lgr(maxint)=1.
         

C     calculate k(g)

      ka(1)=area(1)
      kr(1)=ke1
      ka(ng)=area(maxint)
      kr(ng)=ke2


      do 11 i=2,ng-1
         do 77 j=1,maxint
          if(lgr(j).ge.gr(i))then
           frac=(gr(i)-lgr(j-1))/(lgr(j) - lgr(j-1))
           kr(i)=lkr(j-1)+frac*(lkr(j)-lkr(j-1))
           ka(i)=area(j-1)+0.5*(kr(i) + lkr(j-1))*(gr(i)-lgr(j-1))
           goto 78
          end if
77       continue
78     continue
11    continue


      do 12 i=1,ng-1
       k_g(i)=(ka(i+1)-ka(i))/delg(i)
12    continue

      return

      end
