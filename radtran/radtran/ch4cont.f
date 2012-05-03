************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE CH4CONT
C
C	Returns the 'continuum methane absorption below 1um from
C       Karkoschka (1994), Fink et al above 1um.
C	Wavelengths are supplied in microns.
C
C	Note that since FORTRAN calls pass only the pointers to arrays, 
C	the call to LOCATE supplies only the first nrows elements of
C	TABLE, and since arrays are stored as table(1,1), table(2,1), 
C	etc., then the call to LOCATE passes all the wavelengths only.
C
C-----------------------------------------------------------------------

        real function ch4cont(nu)

	implicit none
	integer		n
	parameter	(n=1810)
	integer		I, J, k, m, k1, k2, nread
	real	lambda, kabs, table(n,2), nu, tmp, frac, xres, xwid
        real  sum1,sum2,wf,dx,x
        character*100 aname
        common /methane/ nread,table

        xres = 0.0125
        xwid = xres/2.0
        if(nread.ne.-999)then
         aname = 'karkoschka.dat'
	 call datarchive(aname)
         open(12,file=aname,status='old')
         print*,'Reading Karkoshka methane data'
         do 100 i=1,n 
          read(12,*)x,table(i,2)
          table(i,1)=0.001*x
100      continue
         nread = -999
         close(12)
        endif

        lambda = 10.0
        if(nu.gt.1000.0)lambda = 1e4/nu

	if ((lambda.lt.table(1,1)).or.(lambda.gt.table(n,1))) then
                ch4cont = 0.0
                return
	endif

	call locate (table, n, lambda, k)

        frac = (lambda - table(k,1))/(table(k+1,1)-table(k,1))
        kabs = table(k,2) + frac*(table(k+1,2)-table(k,2))

        ch4cont = kabs/26850.0	! Convert to units of 1e20 cm-2 

	return

	end
