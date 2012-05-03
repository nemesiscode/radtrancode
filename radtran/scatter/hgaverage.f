      subroutine hgaverage(nlays,ncont,taus,eps,lfrac,tauray,
     1   vwave,f1,g11,g21)
C     ***************************************************************
C     Subroutine to find the average H-G phase function in each layer 
C     in an atmosphere by adding up the OD-weighted HG phase functions
C     of constituent aerosols, and also Rayleigh scattering if required.
C
C     Input variables
C 	nlays	integer	Number of atmospheric layers
C	ncont	integer	Number if aerosol types
C	taus(maxscatlay) real optical depth of each layer
C	eps(maxscatlay) double fraction of thermal emission in each layer
C	lfrac(maxcon,maxscatlay) real fraction of scattering for each 
C				   particle type in each layer
C	tayray(maxscatlay) real optical depth in each layer due to Rayleigh
C			      scattering
C	vwave real wavenumber/wavelength
C
C     Output variables
C	f1(maxscatlay)  real averaged HG f value
C	g11(maxscatlay) real averaged HG g1 value
C	g12(maxscatlay) real averaged HG g2 value
C
C     Pat Irwin		16/11/05
C 
C     ***************************************************************
      implicit none
      include '../includes/arrdef.f'
      integer nlays,ncont,ilay,j,k
      real f1(maxscatlay),g11(maxscatlay),g21(maxscatlay),
     1 taus(maxscatlay)
      real tauray(maxscatlay),lfrac(maxcon,maxscatlay),vwave
      real tf(maxcon),tg1(maxcon),tg2(maxcon)
      double precision f,g1,g2,eps(maxscatlay),omega
      real frac,sfrac,taur,tauscat,tau1

C      print*,'hgaverage: ',nlays,ncont
      do 10 j=1,ncont
       call read_hg(vwave,j,ncont,f,g1,g2)
       tf(j)=sngl(f)
       tg1(j)=sngl(g1)
       tg2(j)=sngl(g2)
10    continue

      do 20 j=1,nlays
       f1(j)=0.0
       g11(j)=0.0
       g21(j)=0.0
C       print*,taus(j),tauray(j),eps(j),(lfrac(k,j),k=1,ncont)
20    continue

      do 100 ilay = 1,nlays
       omega = 1.0 - eps(ilay)
       tauscat = sngl(taus(ilay)*omega)
       taur = tauray(ilay)
C       print*,ilay,omega,taus(ilay),tauscat,taur
C      Calling codes now already include Rayleigh optical depth in
C      tauscat if IRAY=1, so we need to subtract it first here
       tauscat=tauscat-taur

       sfrac=0.0

       do 50 j=1,ncont
        frac = lfrac(j,ilay)
        tau1 = frac*tauscat
        if((tauscat+taur).gt.0.0)then
         frac = tau1/(tauscat+taur)
        else
         frac = 0.0
        endif
        sfrac=sfrac+frac
        if(frac.gt.0)then
         f1(ilay)=f1(ilay)+frac*tf(j)
         g11(ilay)=g11(ilay)+frac*tg1(j)
         g21(ilay)=g21(ilay)+frac*tg2(j)
        endif
50     continue

       if((tauscat+taur).gt.0.0)then
         frac = taur/(tauscat+taur)
       else
         frac = 0.0
       endif

       sfrac=sfrac+frac
C       print*,'frac,sfrac = ',frac,sfrac

       if(frac.gt.0)then
        f1(ilay)=f1(ilay)+frac*0.5
        g11(ilay)=g11(ilay)+frac*0.3
        g21(ilay)=g21(ilay)+frac*(-0.3)
       endif

       if(abs(sfrac-1.0).gt.0.001) then
        print*,'Error in hgaverage sfrac <> 1.0',sfrac
        print*,'Layer = ',ilay
       endif

100   continue

      return

      end
