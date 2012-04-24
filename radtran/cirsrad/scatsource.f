      subroutine scatsource(nlays,nmu,wt,mu,eps,bnu,taus,basehS,
     1 scaleS,f1,g11,g21,nlayerf,basehf,umif,uplf,nf,Jsource)
C     *****************************************************************
C     Subroutine to add plane-parallel calculated scattering to the 
C     source function to simulate scattering observations for a given
C     path 
C
C     Pat Irwin		Original	4/8/05
C			Revised		23/10/07
C     *****************************************************************
      implicit none
c      integer maxlay
      include '../includes/arrdef.f'
      include '../includes/pathcom.f'
      include '../includes/laycom.f'

      double precision eps(maxlay)
      real bnu(maxlay),taus(maxlay),scat
      real umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real radtot,f1(maxlay),g11(maxlay),g21(maxlay)
      real basehf(maxlay),baseh1
      double precision taud,tr,trold
      integer nlays,nmu,ntot,ilay,nf,i,nlayerf,nhalf,j
      real theta0,phi0,sf,f,g1,g2,jsou,dtr,basehS(maxlay)
      real Jsource(maxlay),scaleS(maxlay),basehmin,pi,phi
      parameter(dtr = 3.1415927/180.0)
      double precision mu(maxmu),wt(maxmu)
      logical limbflag

C      print*,nlays,nmu
C      print*,(wt(i),i=1,nmu)
C      print*,(mu(i),i=1,nmu)
C      do i=1,nlays
C       print*,i,eps(i),bnu(i),taus(i),basehS(i),scaleS(i)
C      enddo
C      do i=1,nlays
C       print*,f1(i),g11(i),g21(i)
C      enddo
C      print*,nlayerf
C      do i=1,nlayerf
C       print*,basehf(i)
C      enddo
C      print*,nf


C     Needs major overhaul to deal with sunlight need to specify angles
C     much more carefully and pass in correct geometry to the call to 
C     this routine.

      nhalf = int(0.5*nlays)
      basehmin = 1e10
      pi=3.1415927

      limbflag=.false.
      do i=1,nlays
       if(basehS(i).lt.basehmin)basehmin=basehS(i)
       if(i.lt.nlays)then
        if (basehS(i).eq.basehS(i+1))limbflag=.true.
       endif
      enddo

      print*,'scatsource: limbflag = ',limbflag



C      print*,'umif'
C      do i=1,nlayerf
C        print*,(umif(j,i,1),j=1,nmu)
C      enddo
C      print*,'uplf'
C      do i=1,nlayerf
C        print*,(uplf(j,i,1),j=1,nmu)
C      enddo




      do 20 i=1,nlays

       if(limbflag)then       
        theta0 = asin((RADIUS+basehmin)/(RADIUS+basehS(i)))
        if(i.gt.nhalf)theta0=pi-theta0
       else
        theta0 = acos(1.0/scaleS(i))
       endif
       phi=0.0
C       print*,i,(theta0)*180/3.1415927

       f = f1(i)
       g1 = g11(i)
       g2 = g21(i)

       baseh1 = basehS(i)
 
       if(eps(i).lt.1.0) then
          call fluxintf(nmu,mu,wt,nlayerf,basehf,umif,uplf,
     1    nf,baseh1,f,g1,g2,theta0,phi0,scat)
       else
           scat = 0.0
       endif

C       print*,scat

       Jsource(i)=sngl(eps(i)*bnu(i)+(1.0-eps(i))*scat)
C       print*,i,eps(i),bnu(i),scat,Jsource(i)
20    continue

      return

      end

