      subroutine scatsourcesol(nlays,nmu,wt,mu,eps,bnu,basehS,ibaseH,
     1 tau,emiss_ang,sol_ang,aphi,solar,f1,g11,g21,nlayerf,basehf,umift,
     2 uplft,nf,Jsource,Ig,iwave)
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
      real bnu(maxlay),tau(maxlay),scat
      real umift(maxmu,maxmu,maxscatlay,maxf)
      real uplft(maxmu,maxmu,maxscatlay,maxf)
      real radmij0(181),radmij1(181),radplj0(181),radplj1(181)
      real radmi(maxmu,maxscatlay,181),radpl(maxmu,maxscatlay,181)
      real umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real umip(maxmu,maxscatlay,181),uplp(maxmu,maxscatlay,181)
      real radtot,f1(maxlay),g11(maxlay),g21(maxlay)
      real basehf(maxlay),baseh1,sol_ang,aphi,refZ
      real emiss_ang,xmusun,f0,xf1,xf2,solar
      integer ilayer
      double precision taud,tr,trold
      integer nlays,nmu,ntot,nf,i,nlayerf,nhalf,j,j0,i1,k,iaz
      real theta0,phi0,sf,f,g1,g2,jsou,dtr,basehS(maxlay)
      real tzsat(maxlay),tzsun(maxlay),tphi(maxlay)
      real Jsource(maxlay),basehmin,pi,phi,tausun(maxlay)
      real radsol0,Stheta0,delhS(maxlay)
      integer Ig,iwave,ibaseH(maxlay)

      parameter(dtr = 3.1415927/180.0)
      double precision mu(maxmu),wt(maxmu)
      logical limbflag


      call calc_angle(sol_ang,emiss_ang,aphi,Radius,nlays,baseHS,
     1  ibaseH,nlayer,baseH,delH,tau,tzsat,tzsun,tphi,tausun)

      nhalf = int(0.5*nlays)
      basehmin = 1e10
      pi=3.1415927


      limbflag=.false.
      if(emiss_ang.lt.0.0) limbflag=.true.

      do 20 i=1,nlays

C      Find local solar zenith angle.
C       print*,i,tzsun(i),tzsat(i),tphi(i)
       xmusun = cos(tzsun(i)*dtr) 
       theta0 = tzsat(i)
       Stheta0 = tzsun(i)*dtr
       phi0 = tphi(i)
       radsol0 = solar*exp(-tausun(i))
       j0=-1
C       print*,'i, xmusun = ',i,xmusun
C       print*,'solar,radsol0 = ',solar,radsol0
       if(xmusun.le.sngl(mu(1)))then
        xmusun=sngl(mu(1))
        j0=1
        f0=0.0
       endif
C      Interpolate radiation field to local solar zenith angle.
       do j=1,nmu-1
        if(xmusun.gt.sngl(mu(j)).and.xmusun.le.sngl(mu(j+1)))then
         j0 = j
         f0 = (xmusun-sngl(mu(j)))/sngl((mu(j+1)-mu(j)))
        endif
       enddo
       if(j0.lt.0)then
        print*,'Error in scatsourcesol - j0 not set'
        print*,(mu(j),j=1,nmu),xmusun
        stop
       endif

C       open(12,file='test.dat1',status='unknown')
C        write(12,*)(umift(nmu,nmu,j,1),j=1,nlayerf)
C        write(12,*)(uplft(nmu,nmu,j,1),j=1,nlayerf)
C       close(12)
C       stop

       do j=1,nmu
        do ilayer=1,nlayerf
         do iaz=1,181
          radmij1(iaz)=0.0
          radmij0(iaz)=0.0
          radplj1(iaz)=0.0
          radplj0(iaz)=0.0
         enddo
        
         do i1 = 1,nf+1
          xf1=1.0
          if (i1.gt.1) xf1=2.0
          do iaz=1,181
           xf2 = xf1*cos(float((i1-1)*(iaz-1))*dtr)
           radmij0(iaz)=radmij0(iaz)+umift(j0,j,ilayer,i1)*xf2
           radmij1(iaz)=radmij1(iaz)+umift(j0+1,j,ilayer,i1)*xf2
           radplj0(iaz)=radplj0(iaz)+uplft(j0,j,ilayer,i1)*xf2
           radplj1(iaz)=radplj1(iaz)+uplft(j0+1,j,ilayer,i1)*xf2
          enddo
         enddo

         do iaz=1,181
          radmi(j,ilayer,iaz)=(1.0-f0)*radmij0(iaz)+f0*radmij1(iaz)
          radpl(j,ilayer,iaz)=(1.0-f0)*radplj0(iaz)+f0*radplj1(iaz)
         enddo

        enddo
       enddo
       
       f = f1(i)
       g1 = g11(i)
       g2 = g21(i)


C       if(Ig.eq.5.and.iwave.eq.1)then
C          open(48,file='scatsourcesol.dat',status='unknown',
C     1          form='unformatted')
C                write(48)radmi
C                write(48)radpl
C          close(48)
C       endif

C       stop
       baseh1 = basehS(i)

       if(eps(i).lt.1.0) then
          call fluxintsolA(nmu,mu,wt,nlayerf,basehf,radmi,radpl,
     1    baseh1,f,g1,g2,theta0,Stheta0,phi0,radsol0,scat)
       else
           scat = 0.0
       endif

C       print*,'i,eps,bnu,scat',i,eps(i),bnu(i),scat
       Jsource(i)=sngl(eps(i)*bnu(i)+(1.0-eps(i))*scat)

20    continue
   
C      stop

      return

      end

