      subroutine interp_phase(vv,ncont)
C     ********************************************************************

C	This routine interpolates phase information from the supplied
C	PHASE.DAT files to the current wavenumber/wavelegth. If the
C	variable FIRST is not set to 1, then all the information is read
C	in and stored in a common block for subsequent use. If the arrays
C	are already stored, then Numerical Recipes routine LOCATE is used
C	to find the bracketing values, and linear interpolation is
C	applied.

C     ********************************************************************

	implicit none
	integer		maxpts
	parameter	(maxpts=200)
	integer		ncont, irec1(10), nphas, maxrec(10), I, first,
     1			iunit, npts(10), J, K, npoint, irec, L
	real  		vv, thet(50),phas(10,3,50),head(10,3,3),
     1			cthet(50),inhead(10,maxpts,3),
     2		        inphas(10,maxpts,50),wave(maxpts),
     3                  v0, v1, dv, frac,pi
        parameter (pi=3.1415927)
       	character*512 	buffer
        character*100   ipfile

	common/phase/thet,cthet,phas,head,irec1,nphas,maxrec
	common/interp_phase_initial/first
	common/interpsto/inhead, inphas, npts

 	if (first.ne.1) then
		first = 1
		write (*,*) ' Initialising Interp_phase'
      		do i = 1, ncont
		  iunit = 10+i
                  ipfile='PHASEN.DAT'
                  ipfile(6:6) = char(I+48)
                  open(iunit,file=ipfile,status='old',
     1               access='direct',recl=512, form='formatted')


         	  read(iunit,1000,rec=1)buffer	! Header
         	  if(buffer(2:2).eq.'w')then
            		read(buffer(12:512),*)
     1				v0,v1,dv,npoint,nphas
     			print*,v0,v1,dv,npoint,nphas
		  else
            		read(buffer,*)v0,v1,dv,npoint,nphas
            		print*,v0,v1,dv,npoint,nphas
		  endif

C		  irec = irec1(I)
		  irec=3
	    	  read(iunit,1000,rec=irec) buffer
                  read(buffer,*)(thet(L),L=1,nphas)
                  do L=1,nphas
                   cthet(L)=cos(thet(L)*pi/180.)
                  enddo

		  do J = 1, npoint
			irec = irec + 1
	       		read(iunit,1000,rec=irec) buffer
C                        print*,buffer
C			read(buffer,1010) (inhead(I,J,K),
C     1				K = 1, 3), (inphas(I,J,L),
C     2				L = 1, nphas)

			read(buffer,*) (inhead(I,J,K), 
     1				K = 1, 3), (inphas(I,J,L),
     2				L = 1, nphas)
		  enddo
		  npts(I) = npoint
                  close(iunit)
		enddo
	endif

	do I = 1, ncont
		do J = 1, npts(I)
			wave(J) = inhead(I,J,1)
		enddo
		call locate (wave, npts(I), vv, J)

		if ((J.eq.0).or.(J.eq.npts(I))) then
			iunit = 10 + I
         		print*,'Wavenumber not covered by ',
     1				'scattering file'
         		print*,'Wavenumber = ',vv
         		print*,'Scattering unit number = ',iunit
         		print*,'Permitted range: '
         		read(iunit,1000,rec=1)buffer
         		if(buffer(2:2).eq.'w')then
            			read(buffer(12:512),*)
     1					v0,v1,dv,npoint,nphas
         		else
            			read(buffer,*)v0,v1,dv,npoint,nphas
         		end if
         		print*,v0,v1,dv,npoint,nphas
			stop
		endif

	       	frac=(vv-wave(J))/(wave(J+1)-wave(J))
		head(I,1,1) = wave(J)
		head(I,2,1) = wave(J+1)
	       	head(i,3,1) = vv

		head(I,1,2) = inhead(I,J,2)
		head(I,2,2) = inhead(I,J+1,2)
       		head(i,3,2) = (1.-frac)*head(i,1,2)+frac*head(i,2,2)

		head(I,1,3) = inhead(I,J,3)
		head(I,2,3) = inhead(I,J+1,3)
		head(i,3,3)=(1.-frac)*head(i,1,3)+frac*head(i,2,3)

       		do K=1,nphas
			phas(I,1,K) = inphas(I,J,K)
			phas(I,2,K) = inphas(I,J+1,K)
        	 	phas(i,3,K) = (1.-frac)*phas(i,1,K) +
     1				frac*phas(i,2,K)
		enddo
	enddo

1000	format(a512)
1010    format(1X,F8.2,1X,E12.5,1X,E12.5,1X,50(F10.5,1X))

        return

        end



