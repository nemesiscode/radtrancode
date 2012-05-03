C     $Id: overlap.f,v 1.4 2011-06-17 15:23:00 irwin Exp $
C     ********************************************************************
C     This subroutine combines the absorption coefficient distributions of
C     two overlapping gases. The overlap is implicitly assumed to be random
C     and the k-distributions are assumed to have NG-1 mean values and NG-1
C     weights. Correspondingly there are NG ordinates in total.
C
C     Input variables:
C	delg(ng)	real	Widths of bins in g-space
C	ng		integer Number of ordinates
C	k_g1(ng)	real	k-distribution of first gas
C	q1		real	fraction of first gas
C	k_g2(ng)	real	k-distribution of second gas
C	q2		real	fraction of second gas
C
C     Output variable
C	k_g(ng)		real	k-distribution of both gases together
C
C     Pat Irwin		15/2/95
C
C     ********************************************************************

	subroutine overlap(delg,ng,k_g1,q1,k_g2,q2,k_g)

	implicit none
	integer	ng, I, J, nloop, jtest
	real delg(ng), k_g1(ng),  k_g2(ng), k_g(ng), q1, q2
        include '../includes/arrdef.f'
        real contri(maxrank),weight(maxrank)


C       abort if first k-distribution or abundance = 0.0
        if(k_g1(ng).eq.0.0.or.q1.eq.0.0)then
c         print*,'FIRST CASE q1,q2 = ',q1,q2
c         print*,'k_g1 = ',(k_g1(i),i=1,ng)
         do I = 1, ng
          k_g(I)=k_g2(I)*q2/(q1+q2)
         enddo
         return
        endif

C       abort if second k-distribution or abundance = 0.0
        if(k_g2(ng).eq.0.0.or.q2.eq.0.0)then
c         print*,'SECOND CASE q1,q2 = ',q1,q2
c         print*,'k_g2 = ',(k_g2(i),i=1,ng)
         do I = 1, ng
          k_g(I)=k_g1(I)*q1/(q1+q2)
         enddo
         return
        endif

C       two valid k-distributions read. Need to combine...
        jtest=ng*ng
        if(jtest.gt.maxrank) then 
         print*,'Dimension error in overlap.f'
         print*,'Arrays weight and contri will overflow'
         stop
        endif

	nloop = 0
      	do I = 1, ng
		do J = 1, ng
			nloop = nloop+1
        		weight(nloop) = delg(I)*delg(J)
        		contri(nloop) = (k_g1(I)*q1 + k_g2(J)*q2)/
     1                                  (q1+q2)
 		enddo
	enddo

	call rank(delg, ng, contri, weight, nloop, k_g)

      	return

      	end
