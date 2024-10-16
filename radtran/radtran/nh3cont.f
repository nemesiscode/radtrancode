************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE NH3CONT
C
C	Returns the 'continuum ammonia absorption below 1um from
C       Lutz and Owen (1980).
C
C	Note that since FORTRAN calls pass only the pointers to arrays, 
C	the call to LOCATE supplies only the first nrows elements of
C	TABLE, and since arrays are stored as table(1,1), table(2,1), 
C	etc., then the call to LOCATE passes all the wavelengths only.
C
C-----------------------------------------------------------------------

        real function nh3cont(nu)

	implicit none
	integer		n
	parameter	(n=130)
	integer		I, J, k, m, k1, k2, nread
	real	lambda, kabs, table(n,2), nu, tmp, frac
        real  sum1,sum2,wf,dx,x
        character*100 aname,buffer
        common /ammonia/ nread,table
        integer idiag,iquiet
        common/diagnostic/idiag,iquiet

        if(nread.ne.-999)then
         aname = 'lutzowen.dat'
	 call datarchive(aname)
         open(12,file=aname,status='old')
         if(idiag.gt.0)print*,'Reading Lutz+Owen NH3 data'
         read(12,1)buffer
1        format(a)
         read(12,1)buffer
         read(12,1)buffer
         do 100 i=1,n 
          read(12,*)table(i,1),table(i,2)
100      continue
         nread = -999
         close(12)
        endif

        lambda = 10.0
        if(nu.gt.1000.0)lambda = 1e4/nu

	if ((lambda.lt.table(1,1)).or.(lambda.gt.table(n,1))) then
                nh3cont = 0.0
                return
	endif

	call locate (table, n, lambda, k)

        frac = (lambda - table(k,1))/(table(k+1,1)-table(k,1))
        kabs = table(k,2) + frac*(table(k+1,2)-table(k,2))

C       LutzOwen data in units of m-1 amagat-1
C       1 amagat = 2.6867805e25 molecule/m3
C       Therefore 1 m-amagat = 2.6867805e25 molecule/m2
C       Therefore 1 m-amagat = 2.6867805e21 molecule/cm2 
C       therefore 1 (m-1 amagat-1) = 3.72193e-22 cm2/molecule
C       therefore 1 (m-1 amagat-1) = 3.72193e-2 (1e20) cm2/molecule
C       Check - 1000/26850 = 0.0372439, which is the same.

        nh3cont = 1000*kabs/26850.0	! Convert to units of 1e20 cm2/mol 

	return

	end
