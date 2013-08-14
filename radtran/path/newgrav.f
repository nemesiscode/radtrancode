********************************************************************************
********************************************************************************
C-------------------------------------------------------------------------------
C
C                           SUBROUTINE NEWGRAV
C
C      Sets the gravitational acceleration based on selected planet and 
C      latitude. Gravity is calculated normal to the surface, and corrected 
C      for rotational effects. The input latitude is assumed planetographic.
C
C      The radial vector component of the gravitational acceleration is
C      given by Lindal et al., 1986, Astr. J., 90 (6), 1136-1146, as the
C      expression
C                      ___inf
C          G M  (    \                    (R)**2n                  )
C      g = ---- ( 1 - \     (2n + 1) J_2n (-)     P_2n (sin(latc)) )
C          r**2 (     /__n=1              (r)                      )
C
C                       2
C                    -  - omega**2 r (1 - P_2(sin(latc)))
C                       3
C
C      where r = radial distance, G = gravitational constant, M = planet
C      mass, J_2n is the 2nth zonal harmonic coefficient, P_2n is the
C      Legendre polynomial of order 2n, and latc is the planetocentric
C      latitude. Note that the notation is using spherical coords based
C      on the planet centre. The last part of the expression is simply
C      the centrifugal component and is more simply written as
C
C                           omega**2 r cos(latc)**2.
C
C      The latitudinal vector component is given by the expression
C
C                 ___inf
C          G M  ( \            (R)**2n  d P_2n (sin(latc)) )
C      g = ---- (  \      J_2n (-)      ------------------ )
C          r**2 (  /__n=1      (r)            d latc       )
C
C                      1            d P_2 (sin(latc))
C                    + - omega**2 r -----------------
C                      3                 d latc
C
C      where the notation is as above. The derivatives of the Legendre
C      polynomials can be found through the recurrence relation
C
C                        d P_2n (z)
C             (z**2 - 1) ---------- = 2 n (z P_2n (z) - P_(2n-1) (z)).
C                           d z
C
C      Using this relation, it can be seen that the last expression of
C      the latitudinal component represents centrifugal acceleration and
C      is more clearly written as
C
C                           omega**2 r cos(latc) sin(latc).
C
C      These two gravitational vector components are normal to each other
C      and must be summed to give the overall gravity acting normal to the
C      planetary surface.
C
C      Original version    A.Weir
C      Adapted for general use 29/4/96  Pat Irwin
C-------------------------------------------------------------------------------

       subroutine newgrav (iplanet, lat_in, h, radius, g, pname) 
C      *******************************************************************
C      Input variables:
C	iplanet		integer		Planet ID code
C	lat		real		Planetocentric Latitude (degrees)
C	h		real		Height above reference surface (km)
C
C      Output variables:
C       radius		real		Planetary radius (km) at latitude lat
C	g		real		gravity acc. (m/s2)
C	pname		character*8	Planet name
C
C      *******************************************************************

       implicit none
       INCLUDE '../includes/planrad.f'
       integer       iplanet, I, isurf
       real        lat, g, pi, latc, thetagc, slatc, r, legpol,
     2               clatc, gtheta, pol(6), Rr, radius,lat_in
       real xgm,xcoeff(3),xradius,xellip,xomega,h
       character*8 pname
       data pi/3.141592654/ 

C-------------------------------------------------------------------------------
C
C      First read in the gravitometric data
C
C-------------------------------------------------------------------------------

       call read_grav(iplanet,xgm,xcoeff,xradius,xellip,xomega,pname,
     1   isurf)

C-------------------------------------------------------------------------------
C
C      First convert to planetocentric latitude. Then calculate radial
C      distance to the centre of mass. Precalculate some things for speed.
C
C-------------------------------------------------------------------------------

       lat = 2 * pi * lat_in/360.
       latc = thetagc(lat, xellip)
       slatc = sin(latc)
       clatc = cos(latc)
       
       if(jradf.gt.0)xradius=radius2*1.e5

C      Rr is the ratio of radius at equator to radius at current latitude
       Rr = sqrt(clatc**2 + (xellip**2 * slatc**2))
       r = (xradius+h*1e5)/Rr
C       radius = xradius*1e-5
       radius = (xradius/Rr)*1e-5

       do I = 1, 6
              pol(I) = legpol(I,slatc)
       enddo

C-------------------------------------------------------------------------------
C
C      Evaluate radial contribution from summation for first three terms,
C      then subtract centrifugal effect.
C
C-------------------------------------------------------------------------------

       g = 1.
       do I = 1, 3
              g = g - ((2*I+1) * Rr**(2 * I) * xcoeff(I) * pol(2*I))
C              print*,g*0.01*xgm/r**2
       enddo

       g = (g * xgm/r**2) - (r * xomega**2 * clatc**2)
C       print*,'gradial = ',g*0.01
C-------------------------------------------------------------------------------
C
C      Evaluate latitudinal contribution for first three terms, then add
C      centrifugal effects.
C
C-------------------------------------------------------------------------------

       gtheta = 0.
       do I = 1, 3
              gtheta = gtheta - (4 * I**2 * Rr**(2 * I) * xcoeff(I) 
     1               * (pol(2*I-1) - slatc * pol(2*I))/clatc)
C       print*,gtheta*0.01*xgm/r**2
       enddo

       gtheta = (gtheta * xgm/r**2) + (r * xomega**2 * 
     1        clatc * slatc)

C       print*,'gtheta = ',gtheta*0.01
C-------------------------------------------------------------------------------
C
C      Combine the two components and write the result.
C
C-------------------------------------------------------------------------------

       g = sqrt(g**2 + gtheta**2)*0.01
C       print*,'gravity =',g
       return       

       end      

********************************************************************************
********************************************************************************

C-------------------------------------------------------------------------------
C
C                           FUNCTION THETAGC
C
C      Converts planetographic latitude to planetocentric. Angles must be
C      supplied in radians.
C
C-------------------------------------------------------------------------------

       function thetagc (lat, e)

       implicit none
       real        e, lat, thetagc

       thetagc = atan(tan(lat)/e**2)

       end

********************************************************************************
********************************************************************************
C-------------------------------------------------------------------------------
C
C                           FUNCTION LEGPOL
C
C      Calculates zero order Legendre polynomials based on the recurrence
C      relation 
C               (n) P   (z) = (2n-1) z P    (z) - (n-1) P     (z).
C                    (n)                (n-1)            (n-2)
C
C-------------------------------------------------------------------------------

       function legpol (n, z)
       
       implicit none
       integer       I, n
       real        legpol, z, pol(3)

       pol(1) = 1.
       pol(2) = z

       if ((n.eq.0).or.(n.eq.1)) then
              legpol = pol(n+1)
       else
              I = 1
10            I = I + 1

              if (I.gt.50) then
                     write (*,*) '* * * * * * * * * * * * *'
                     write (*,*) '*                       *'
                     write (*,*) '*   Looping in LEGPOL   *'
                     write (*,*) '*   Program terminated  *'
                     write (*,*) '*                       *'
                     write (*,*) '* * * * * * * * * * * * *'
              endif

              pol(3) = (((2 * I - 1) * pol(2) * z) 
     1               - ((I - 1) * pol(1)))/float(I)
              if (I.eq.n) then
                     legpol = pol(3)
                     return
              else
                     pol(1) = pol(2)
                     pol(2) = pol(3)
                     goto 10
              endif
       endif


       end

********************************************************************************
********************************************************************************
