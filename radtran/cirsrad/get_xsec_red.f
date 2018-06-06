************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				 SUBROUTINE GET_XSEC
C
C	Rearranges variables for use in common block dust as opposed to 
C	dustspec. 
C
C-----------------------------------------------------------------------

	subroutine get_xsec_red(xscfil,ncont,vwave)

	implicit none

	include '../includes/arrdef.f'

	integer		ncont, nsec, ncont1, nsec1, I, J, icont
	real		xsec(2,maxcon,maxsec), vsec(maxsec), 
     1   vsec1(maxsec), xsec1(2,maxcon,maxsec)
	real            vwave
	double precision f, g1, g2, asym
	character*100	xscfil

	common/dustspec/xsec,vsec,nsec
	common/dust/vsec1,xsec1,nsec1,ncont1

C-----------------------------------------------------------------------
C
C	Read from XSC file, then assign values.
C
C-----------------------------------------------------------------------

	call read_xsec(xscfil)

        print*,'get_xsec_red: read_xsec OK'
	ncont1 = ncont
	nsec1 = nsec
        print*,'nsec1 = ',nsec
	do I = 1, nsec
		vsec1(I) = vsec(I)
		do J = 1, ncont
			call read_hg(vwave,J,ncont,f,g1,g2)			
                        print*,'get_xsec_red ',f,g1,g2
			asym = f*g1+(1-f)*g2
			print*,'get_xsec_red ',f, g1, g2
     			print*,'get_xsec_red ',asym
			if (asym.GE.0) then
			  xsec1(1,J,I) 
     1	                  = xsec(1,J,I)-real(asym)*xsec(2,J,I)
		          print*,'get_xsec ',real(asym)
			  print*,'get_xsec ',xsec(2,J,I)
		          print*,'get_xsec ',real(asym)*xsec(2,J,I)
	                  xsec1(2,J,I) = xsec(2,J,I)
			else
				xsec1(1,J,I) = xsec(1,J,I)
				xsec1(2,J,I) = xsec(2,J,I)
			endif
			print*,'get_xsec_red xsec1, xsec'
                        print*,xsec1(1,J,I),xsec(1,J,I)
		enddo
	enddo

C-----------------------------------------------------------------------
C
C	Return and end.
C
C-----------------------------------------------------------------------

	return

	end
			
************************************************************************
************************************************************************
