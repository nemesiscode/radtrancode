************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				SUBROUTINE RANK 
C
C     	Pat Irwin		23/2/96
C	A.L.Weir	02/03/96
C
C-----------------------------------------------------------------------

	subroutine rank(gw,ng,cont,weight,nloop,k_g)
	implicit none
        include '../includes/arrdef.f'

	integer		ng, nloop, I, ig
	real		gw(ng),cont(nloop), weight(nloop), k_g(maxg), 
     1			g_ord(maxg), gdist(0:maxrank), sum, frac

C-----------------------------------------------------------------------
C
C	Sum delta gs to get cumulative g ordinate. rank cont and weight
C	arrays in ascending order of k (i.e. cont values).
C
C-----------------------------------------------------------------------

C        print*,'cirsrad rank'
	g_ord(1)=0.
	do I = 1, ng
		g_ord(I+1) = g_ord(I) + gw(I)
	end do

	call sort2(nloop, cont, weight)

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
	enddo

	ig = 1
	sum=0.
	do I = 1, nloop
		if(gdist(I).lt.g_ord(ig+1))then
			k_g(ig) = k_g(ig) + cont(I) * weight(I)
			sum = sum + weight(I)
		else
			frac = (g_ord(ig+1)-gdist(I-1))/
     1				(gdist(I)-gdist(I-1))
			k_g(ig) = k_g(ig) + frac * cont(I)*weight(I)
			sum = sum + frac * weight(I)
			k_g(ig) = k_g(ig)/sum
			ig = ig + 1
			sum = (1.-frac) * weight(I)
			k_g(ig) = k_g(ig) + (1.-frac)*cont(I)
     1				*weight(I)
		endif
	enddo
	if (ig.eq.ng) k_g(ig) = k_g(ig)/sum

	return

	end

************************************************************************
************************************************************************
