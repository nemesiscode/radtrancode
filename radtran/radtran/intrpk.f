************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				SUBROUTINE INTRPK
C
C	Gets estimate of g from k using linear interpolation. Called
C	by KML. Note that interpolation is done in multiplicative steps
C	rather than log steps, but that logs are not explicitly used 
C	until the last moment to speed things up.
C
C	Bracketing values of k are subdivided into intervals, and
C	g evaluated at each point. They are then searched and a new
C	k value interpolated. The new bracketing values are returned in 
C	place of the old.	
C
C-----------------------------------------------------------------------

	subroutine intrpk(kb,kt,aa,bb,BL,n,gin,kout)

	integer		max
	parameter	(max=1000)
	integer	nint, I
	double precision kb, kt, gin, kout, dlk, lkr(max), k1, 
     1    lgr(max), slvf, frac, lkt, lkb, gb, gt, aa, bb, BL

C-----------------------------------------------------------------------
C
C	Check inputs are OK
C
C-----------------------------------------------------------------------

	if (n.gt.max) then
		write (*,*) ' Problems in INTRPK '
		write (*,*) ' Nint exceeds max: ', n, max
		stop
	endif

	gb = dslvf(kb,aa,bb,BL,0.d0)
	gt = dslvf(kt,aa,bb,BL,0.d0)
	if ((gin.gt.gt).or.(gin.lt.gb)) then
	   write (*,*)' INTRPK: Supplied Ks do not bracket target: '
	   write (*,*) ' Target: ', gin
	   write (*,*) ' k: ', kb, ' gives ',gb
	   write (*,*) ' k: ', kt, ' gives ',gt
	   stop
	endif
	if (gb.eq.gin) then 
		kout = kb
		goto 10
	elseif (gt.eq.gin) then
		kout = kt
		goto 10
	endif
		
C-----------------------------------------------------------------------
C
C	Do interpolation
C
C-----------------------------------------------------------------------

	dlk = (kt/kb)**(1./real(n-1))	
	lkr(1)=kb
        do I = 2, n
                k1 = lkr(I-1)*dlk
                lkr(I) = k1
                lgr(I) = dslvf(k1,aa,bb,BL,0.d0)
        enddo

        if (gin.eq.0.) then
                kout=kb
        else if (gin.eq.1.) then
                kout=kt
        else
                do I = 1, n
                        if (lgr(I).ge.gin) then
                                frac = (gin-lgr(I-1))/(lgr(I)-lgr(I-1))
                                kt = lkr(I)
                                kb = lkr(I-1)
                                lkt = dlog(kt)
                                lkb = dlog(kb)
                                kout = lkb + frac*(lkt-lkb)
                                kout = dexp(kout)
                                goto 10
                        endif
                enddo
10      endif

	return

	end

************************************************************************
************************************************************************


