************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				SUBROUTINE KML
C
C 	For a given input line strength (SL) and line broadening (BL)
C	parameter, returns the value of k at the supplied value of g.
C
C	LO refers to Lacis and Oinas, A description of the Correlated k
C	distribution for modelling nongray gaseous absorption, thermal
C	emission, and multiple scattering in vertically inhomogeneous 
C	atmospheres. JGR 96, No. D5, pp9027-9063, May 1991.	
C
C-----------------------------------------------------------------------

      	subroutine kml(SL,BL,gin,kout)

      	implicit none
	integer		maxint
	real		step
	parameter	(maxint=200,step=0.01)
	integer	count, I
	real	gin, kout, kmax, ke1, ke2, gg, slvf, dlk, k1, 
     1		lkr(maxint), lgr(maxint), frac, 
     2		gout, dgdk
	double precision SL, BL, aa, bb, pi, dg1, dg2, dslvf, kb,
     1    kt, dgg, dkout, dgout, dke1, dke2, tmp, x

	pi = 2.0d0 * dasin(1.0d0)

C-----------------------------------------------------------------------
C
C	First check for trivial case of k = 0. Occurs for SL or BL=0,
C	or if max value = 0. Then define some parameters as used 
C	by LO p9037. Calculate the maximum value of k in the frequency
C	distribution (L&O eq. 19)
C
C-----------------------------------------------------------------------

	if ((SL.eq.0.).or.(BL.eq.0.)) then
		print*,'Kml: Error 0',SL,BL
       		kout=0.
       		return
	endif

      	aa = 0.5*dsqrt(pi*BL*SL)
      	bb = 0.5*dsqrt(pi*BL/SL)
    
       
        x = (pi*BL/3.)**2
        if(abs(x).lt.1e-10)then
         tmp = 0.5*x - 0.125*x**2
        else
         tmp = (dsqrt((pi*BL/3.)**2 + 1.) - 1.)
        endif
      	tmp = (3.*SL/(pi*BL))*tmp
        kmax = sngl(tmp)
      	if (kmax.eq.0.) then
                print*,'kml: error 1',tmp,kmax
       		kout=0.
       		return
      	end if

C-----------------------------------------------------------------------
C
C     	First find rough bounds on the range of k. Starting from kmax, 
C	halve ke1 until g(ke1) - gin is negative. Similarly find ke2
C	s.t. g(ke2) - gin is positive. WE have now bracketed the target.
C
C	Note that g(k) - gin is supplied by function SLVF
C
C-----------------------------------------------------------------------

c      	ke1 = kmax
c      	ke2 = kmax
c      	gg = gin

      	dke1 = dble(kmax)
      	dke2 = dble(kmax)
       	dgg = dble(gin)

	count = 0
c10    	if (slvf(ke1,aa,bb,BL,gg).gt.0) then
10    	if (dslvf(dke1,aa,bb,BL,dgg).gt.0.0d0) then
c      		ke1 = ke1*0.5
       		dke1 = dke1*0.5d0
       		count = count + 1
       		if (count.lt.500) then
        		goto 10
       		else
        		write (*,*) 'Error. ke1 has been halved more ',
     1				'than 500 times'
        		write (*,*) 'dke1,aa,bb,BL,dgg'
        		write (*,*) dke1,aa,bb,BL,dgg
        		stop
       		end if
      	end if

	count = 0
c20    	if (slvf(ke2,aa,bb,BL,gg).lt.0) then
20    	if (dslvf(dke2,aa,bb,BL,dgg).lt.0) then
c        	ke2 = ke2*2.
        	dke2 = dke2*2.0d0
        	count = count+1 
        	if (count.lt.500) then
         		goto 20
        	else
         		write (*,*) 'Error. ke2 has been doubled more ',
     1				'than 500 times'
         		write (*,*) 'dke2,aa,bb,BL,dgg'
         		write (*,*) dke2,aa,bb,BL,dgg
         		stop
       		end if
      	end if

C-----------------------------------------------------------------------
C
C     	Split up range ke1 to ke2 into maxint intervals and determine 
C     	(maxint+1) ordinates and determine for each ordinate k, g(k). 
C	Values go into lkr(i) and lgr(i) Note that intervals are
C	multiplicative raher than additive. Force the extremes of lgr
C	to 0. and 1.
C
C	Do linear interpolation (on logs) to get estimate of kout. Note
C	that we do not use log(k) before actual interpolation to speed
C	things up. 
C
C	By using gg = 0., SLVF now returns g(k).
C
C-----------------------------------------------------------------------

c	kb = dble(ke1)
c	kt = dble(ke2)
c	dgg = dble(gin)
	kb = dke1
	kt = dke2
	call intrpk(kb,kt,aa,bb,BL,maxint,dgg,dkout)
	kout = sngl(dkout)


C-----------------------------------------------------------------------
C
C	Now have a good estimate of k. Use Newtonian iteration to
C	improve it. In pathological cases (usually very large values of
C	BL where the gk curve is steep - see LO fig 5b) this can fail.
C	In these cases, we brute force down to an acceptable value using
C	the above interpolation method.
C
C-----------------------------------------------------------------------

      	gg=gin
	dgg = dble(gin)
	count = 0
      	gout = slvf(kout,aa,bb,BL,gg)

40	ke1 = kout - step*kout
	ke2 = kout + step*kout
   	dke1 = dble(ke1)
      	dke2 = dble(ke2)
      	dg1 = dslvf(dke1,aa,bb,BL,dgg)    
      	dg2 = dslvf(dke2,aa,bb,BL,dgg)
      	dgdk = sngl(dg2 - dg1)/(2.*step*kout)
	if (dgdk.eq.0.) then				! Brute force
50     		count = count+1
		call intrpk(kb,kt,aa,bb,BL,maxint,dgg,dkout)
      		dgout = dslvf(dkout,aa,bb,BL,dgg)
		if (dgout.lt.0.) then
			kb = dkout
		else
			kt = dkout
		endif
      		if ((abs(dgout).gt.1e-4).and.(count.lt.500)) goto 50
	else						! Newtonian	
      		count = count+1
  		kout = kout - gout/dgdk
      		if (kout.lt.ke1) kout = ke1		! Hard stops
      		if (kout.gt.ke2) kout = ke2
      		gout = slvf(kout,aa,bb,BL,gg)
      		if ((abs(gout).gt.1e-4).and.(count.lt.500)) goto 40
	endif

C-----------------------------------------------------------------------
C
C	Return and end
C
C-----------------------------------------------------------------------

      return

      end

************************************************************************
************************************************************************
