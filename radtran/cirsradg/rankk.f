************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				SUBROUTINE RANKK 
C
C     	Pat Irwin	23/2/96
C	A.L.Weir	02/03/96
C	n teanby	7/8/03	copied from cirsrad/rank. oxcirsg uses this 
C				routine and i've renamed it to avoid confusion 
C				with other rank.f's (like the numerical 
C				recipies one). 
C	Pat Irwin	29/2/12	Updated for Radtrans2.0
C
C-----------------------------------------------------------------------

	subroutine rankk(gw,ng,cont,dcont,weight,nloop,k_g,dkdt)

	implicit none
        include '../includes/arrdef.f'

	integer		ng, nloop, I, ig,ico(MAXRANK)
	real		gw(ng),cont(nloop), weight(nloop), k_g(MAXG), 
     1			g_ord(MAXG), gdist(0:2*MAXG), sum, frac,
     2			dkdt(MAXG),dcont(nloop),tmp(MAXRANK)

C-----------------------------------------------------------------------
C
C	Sum delta gs to get cumulative g ordinate. rank cont and weight
C	arrays in ascending order of k (i.e. cont values).
C
C-----------------------------------------------------------------------

	g_ord(1)=0.
	do I = 1, ng
		g_ord(I+1) = g_ord(I) + gw(I)
	end do

        do I=1,nloop
         ico(I) = I
        enddo

C Sort random k-coeffs into order. Integer array ico records which swaps
C have been made so that we can also re-order the gradients.
        CALL sort2g(nloop, cont, ico)

C Resort the weights:
        DO I=1,nloop
          tmp(I) = weight(ico(I))
        ENDDO
        DO I=1,nloop
          weight(I) = tmp(I)
        ENDDO

C Resort the gradient w.r.t. Temperature
        DO I=1,nloop 
          tmp(I) = dcont(ico(I))
        ENDDO
        DO I=1,nloop
          dcont(I) = tmp(I)
        ENDDO

C-----------------------------------------------------------------------
C
C	Now form new g(k) by summing over weight. The new k(g) is then
C	found by ascending the arranged ks and getting the weighted 
C	averages of the k values in the relevant g interval. Linear
C	interpolation is used to deal with steps crossing over the 
C	required g intervals.
C
C-----------------------------------------------------------------------

        gdist(0) = 0.
	gdist(1) = weight(1)
	do I = 2, nloop
		gdist(I) = weight(I) + gdist(I-1) 
	enddo

	do I = 1, ng
		k_g(I) = 0.
		dkdt(I) = 0.
	enddo

	ig = 1
	sum=0.
	do I = 1, nloop
		if(gdist(I).lt.g_ord(ig+1))then
			k_g(ig) = k_g(ig) + cont(I) * weight(I)
			dkdt(ig) = dkdt(ig) + dcont(I)*weight(I) 
			sum = sum + weight(I)
		else
			frac = (g_ord(ig+1)-gdist(I-1))/
     1				(gdist(I)-gdist(I-1))
			k_g(ig) = k_g(ig) + frac * cont(I)*weight(I)
			dkdt(ig) = dkdt(ig) + frac * dcont(I)*weight(I)
			sum = sum + frac * weight(I)
			k_g(ig) = k_g(ig)/sum
			dkdt(ig) = dkdt(ig)/sum
			ig = ig + 1
			sum = (1.-frac) * weight(I)
			k_g(ig) = (1.-frac)*cont(I)*weight(I)
			dkdt(ig) = (1.-frac)*dcont(I)*weight(I)
		endif
	enddo
	if (ig.eq.ng) then
          k_g(ig) = k_g(ig)/sum
          dkdt(ig) = dkdt(ig)/sum
        endif

	return

	end

************************************************************************
************************************************************************
