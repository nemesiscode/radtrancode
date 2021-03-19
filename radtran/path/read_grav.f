      subroutine read_grav(iplanet,xgm,xcoeff,xradius,
     & xellip,xomega,pname1,isurf)
C     $Id: read_grav.f,v 1.2 2004-11-25 11:53:02 irwin Exp $
C     **************************************************************
C     Subroutine to read in planetary gravity data from 'gravity.dat'
C
C     Data is extracted from Astronomical Almanac 1994, E88
C
C     Pat Irwin		Documented	13/4/00
C
C     **************************************************************
      implicit none
      integer mplanet,nplanet,iplanet,isurf,i,j
      real xgm,xcoeff(3),xradius,xellip,xomega,pi,Grav
      parameter (mplanet=87,pi=3.1415927,Grav=6.672E-11)
      real plan_mass(mplanet),Jcoeff(3,mplanet),aradius(mplanet)
      real flatten,rotation,ellip(mplanet),omega(mplanet)
      integer isurftab(mplanet)
      character*100 xbuffer,buffer,aname
      character*8 pname(mplanet),pname1
      include '../includes/planrad.f'
      common /gravdata/nplanet,pname,plan_mass,Jcoeff,aradius,
     1 isurftab,ellip,omega
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


1     format(a)

C      print*,'pname(1) = ',pname(1)
      if(pname(1).ne.'Mercury ')then      
       aname = 'gravity.dat'
       call datarchive(aname)
       if(idiag.gt.0)print*,'Reading gravity data from ',aname
C        print*,'Reading gravity data from ',aname

       open(12,file=aname,status='old')
       read(12,*)nplanet
       if(nplanet.gt.mplanet)then
        print*,'Error in read_grav.f. nplanet > mplanet'
        print*,nplanet,mplanet
        stop
       endif
       do 10 i=1,4
        read(12,1)xbuffer
10     continue
       do 20 i=1,nplanet
        read(12,1)xbuffer
C       Ignore first 4 spaces which are reserved for planet number
        pname(i)=xbuffer(5:12)
        buffer=''
        buffer(1:72)=xbuffer(13:84)
        read(buffer,*)plan_mass(i),(Jcoeff(j,i),j=1,3),aradius(i),
     &   flatten,rotation,isurftab(i)
        ellip(i)=1.0/(1.0-flatten)
        omega(i) = 2*pi/(rotation*24.0*3600.0)
        plan_mass(i)=plan_mass(i)*Grav*1e24*1e6
        Jcoeff(1,i)=Jcoeff(1,i)/1e3
        Jcoeff(2,i)=Jcoeff(2,i)/1e6
        Jcoeff(3,i)=Jcoeff(3,i)/1e8
C       Convert radius to units of cm
        aradius(i)=aradius(i)*1e5
20     continue
       close(12)


       if(iplanet.gt.nplanet)then
        print*,'iplanet is greater than the number of planets listed'
        print*,'in the gravity.dat file:'
        print*,'iplanet, nplanet = ',iplanet,nplanet
        print*,'Check planet ID or add data to gravity.dat'
        print*,'Aborting...'
        stop
       endif

      endif

      xgm = plan_mass(iplanet)
      if(MCMCflag.eq.1)then
       xgm=MCMCmass*Grav*1e24*1e6
      endif
      if(mass2.gt.0.0.and.jloggf.gt.0)then
        xgm = mass2*Grav*1e24*1e6
      else
        mass2 = plan_mass(iplanet)/(Grav*1e24*1e6)
      endif


      do 30 j=1,3
       xcoeff(j)=Jcoeff(j,iplanet)
30    continue
      xradius=aradius(iplanet)
      if(MCMCflag.eq.1)then
       xradius=MCMCrad*1e5
      endif
      if(radius2.gt.0.0.and.jradf.gt.0)then
       xradius=radius2*1e5
      else
       radius2=xradius/1e5
      endif
      
      xellip = ellip(iplanet)
      xomega = omega(iplanet)
      pname1 = pname(iplanet)
      isurf = isurftab(iplanet)
      
      return
      end
