      program readphase
C     ********************************************************************     
C 
C	Program to read in a PHASE*.DAT file

C     ********************************************************************   

	implicit none
	integer		maxpts
	parameter	(maxpts=200)
	integer		ncont, irec1(10), nphas, maxrec(10), I, first, 
     1			iunit, npts(10), J, K, npoint, irec, L 
	real vv, thet(50),phas(10,3,50),head(10,3,3),cthet(50),
     1 inhead(10,maxpts,3), inphas(10,maxpts,50), 
     2 wave(maxpts), v0, v1, dv, frac
      	character*512 	buffer
        character*100 opfile


        call prompt('Enter file name : ')
        read(5,1)opfile
1       format(a)

        print*,'opening: ',opfile
        open(11,file=opfile,status='old',access='direct',
     1    recl=512,form='formatted')
 
        iunit=11

        read(iunit,1000,rec=1)buffer	! Header
        if(buffer(2:2).eq.'w')then
            read(buffer(12:512),*) v0,v1,dv,npoint,nphas
        else
            read(buffer,*)v0,v1,dv,npoint,nphas
	endif
        print*,v0,v1,dv,npoint,nphas

        maxrec(1)=npoint+3
        irec1(1)=3
        irec=3
        read(iunit,1000,rec=irec)buffer
        read(buffer,*)(thet(j),j=1,nphas)
        do j=1,nphas
         cthet(j)=cos(thet(j)*0.0174532)
         print*,thet(j),cthet(j)
        end do

        I=1

        call prompt('Enter required wavenumber : ')
        read*,vv

	irec = irec1(I)
	do J = 1, npoint
	 irec = irec + 1
       	 read(iunit,1000,rec=irec) buffer
	 read(buffer,1010) (inhead(I,J,K), 
     1			    K = 1, 3), (inphas(I,J,L), 
     2			    L = 1, nphas)
	enddo
	npts(I) = npoint	

	do J = 1, npts(I)
	 wave(J) = inhead(I,J,1)
	enddo

	call locate (wave, npts(I), vv, J)
 
	if ((J.eq.0).or.(J.eq.npts(I))) then
	 iunit = 10 + I		
         print*,'Wavenumber not covered by scattering file'
         print*,'Wavenumber = ',vv
         print*,'Scattering unit number = ',iunit 
         print*,'Permitted range: '
         read(iunit,1000,rec=1)buffer
         if(buffer(2:2).eq.'w')then
            read(buffer(12:512),*)v0,v1,dv,npoint,nphas
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


        close(iunit)
        
        call prompt('Enter output filename : ')
        read(5,1)opfile
        open(12,file=opfile,status='unknown')
         write(12,*)'Wavenumber = ',vv
         do k=1,3
          write(12,*)k,(head(I,J,k),J=1,3)
         end do
         write(12,*)nphas
         do k=1,nphas
          write(12,*)thet(k),(phas(I,J,K),J=1,3)
         enddo
        close(12)


1000	format(a512)
1010    format(1X,F8.2,1X,E12.5,1X,E12.5,1X,50(F10.5,1X))


        end


