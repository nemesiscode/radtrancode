      subroutine limbscatter(nlays,nmu,wt,mu,eps,bnu,taus,
     1 f1,g11,g21,umif,uplf,nf,radtot)
C     *****************************************************************
C     Subroutine to add plane-parallel calculated scattering to the 
C     source function to simulate limb scattering observations.
C
C     Pat Irwin		Original	4/8/05
C
C     *****************************************************************
      implicit none
c      integer maxlay
      include '../includes/arrdef.f'
      include '../includes/pathcom.f'
      include '../includes/laycom.f'

      real bnu(maxscatlay),taus(maxlay),scat
      real umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real radtot,f1(maxlay),g11(maxlay),g21(maxlay)
      double precision taud,tr,trold,eps(maxscatlay)
      integer nlays,nmu,ntot,jlay,ilay,nf,i
      real theta0,phi0,sf,f,g1,g2,jsou,dtr
      parameter(dtr = 3.1415927/180.0)
      double precision mu(maxmu),wt(maxmu)
      
      ntot = 1+nlays-botlay
      radtot=0.0
      trold=1.0
      taud=0.0

      do 100 ilay = 1,2*ntot
        if(ilay.le.ntot)then
         jlay = ilay
        else
         jlay = 1+2*ntot-ilay
        endif


        call calczenscale(ilay,nlays,ntot,theta0,
     1   phi0,sf)

        taud = taud + taus(jlay)*sf
        tr = dexp(-taud)

        f = f1(jlay)
        g1 = g11(jlay)
        g2 = g21(jlay)

        if(eps(jlay).lt.1.0) then
          call fluxint(nmu,mu,wt,umif,uplf,nf,jlay,f,g1,g2,
     1     theta0,phi0,scat)
        else
           scat = 0.0
        endif

        jsou =  sngl(eps(jlay)*bnu(jlay) + (1.0-eps(jlay))*scat)

        radtot = radtot + jsou*(trold-tr)
        trold = tr

100   continue

      return

      end

      subroutine calczenscale(ilay,nlays,ntot,
     1  theta0,phi0,sf)
C     *****************************************************************
C     Subroutine to calculate appropriate layer scaling factors for a
C     limb path and also the viewing zenith angle.
C
C     Pat Irwin		Original	4/8/05
C
C     *****************************************************************
      implicit none
 
      include '../includes/arrdef.f'
      include '../includes/pathcom.f'
      include '../includes/laycom.f'

      real theta0,phi0,sf,pi,sin2a,cosa,z0,s1,s0,stmp
      integer ilay,ntot,jlay,jlay1,nlays
      parameter(pi=3.1415927)

      if(ilay.le.ntot)then
         jlay1 = ilay
      else
         jlay1 = 1+2*ntot-ilay
      endif

      jlay = nlays+1-jlay1

      z0 = radius + baseh(botlay)
      sin2a = 1.0
      cosa = 0.0

      stmp= (radius+baseh(jlay))**2 - sin2a*z0**2
C      Sometimes, there are rounding errors here that cause the program
C      to try to take the square root of a _very_ small negative number.
C      This quietly fixes that.
      if(stmp.lt.0.0) stmp=0.0
      s0 = sqrt(stmp) - z0*cosa
      if(jlay.eq.nlays)then
        s1=sqrt((radius+h(npro))**2 - sin2a*z0**2) - z0*cosa
        sf = (s1-s0)/(h(npro)-baseh(jlay))
      else
        s1 = sqrt((radius+baseh(jlay+1))**2 - sin2a*z0**2) - z0*cosa 
        sf = (s1-s0)/(baseh(jlay+1)-baseh(jlay))
      endif

      theta0 = asin((radius+baseh(botlay))/(radius+baseh(jlay)))
      
      if(ilay.gt.ntot)theta0 = pi-theta0
      phi0 = 0.0

      return
       
      end
