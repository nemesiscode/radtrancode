c-----------------------------------------------------------------------
	program interp_prf
c-----------------------------------------------------------------------
c	logarithmic interpolation of a radtran style .prf or .ref file
c-----------------------------------------------------------------------
c	N. Teanby	19-04-04	Original code
c-----------------------------------------------------------------------
      
      implicit none

      integer nvmr_max,npro_max,nvmr,npro,npro1
      parameter (nvmr_max=30,npro_max=300)
      integer iform,iplanet,idiso(nvmr_max,2)
      real latitude,molwt
      real vmrs(npro_max,nvmr_max),vmrs1(npro_max,nvmr_max)
      real vmrs_in(npro_max),vmrs_out(npro_max)
      real height(npro_max),height1(npro_max)
      real press(npro_max),press1(npro_max)
      real press_log(npro_max),press1_log(npro_max)
      real temp(npro_max),temp1(npro_max)
      real p1,p2
      character*100 prffile,prffile_out
      character*1 ans
      integer i,j,l

c  ** get input filename **      
	print*,'Enter input file name (.prf format):'
      read*,prffile        

c  ** read .prf file **    
	call zreadprf(prffile,nvmr_max,npro_max,
     $  iplanet,latitude,nvmr,npro,molwt,idiso,height,press,temp,vmrs)
	print*,'npro = ',npro
	print*,'p1 = ',press(1)
	print*,'p2 = ',press(npro)
c  ** get new number of points
	print*,'Enter new number of points in profile (npro):'
      read*,npro1
	print*,'change pressure range? [y/n]'
      read*,ans
      if ((ans.eq.'Y').or.(ans.eq.'y')) then
        print*,'Enter new bottom and top pressures (p1 and p2):'
        read*,p1,p2
	else
        p1=press(1)
        p2=press(npro)
      endif

c  ** compute log pressure **
	do i=1,npro
        press_log(i) = log(press(i))
      enddo
      
c  ** create new pressure baseline, linearly spaced in log(p) **
      do i=1,npro1
        press1_log(i) = log(p1) + (log(p2)-log(p1))*(i-1)/(npro1-1)
        press1(i)     = exp( press1_log(i) )
      enddo

c  ** linear interp height in log pressure **
	call zlinint(npro,npro_max,press_log,height,npro1,npro_max,
     >  press1_log,height1)

c  ** log interp of temperature **
	call zlinint(npro,npro_max,press_log,temp,npro1,npro_max,
     >  press1_log,temp1)
      
c  ** linearly interpolate each gas profile in log pressure **
	do i=1,nvmr
        do j=1,npro
          vmrs_in(j)=vmrs(j,i)
	  enddo
	  call zlinint(npro,npro_max,press_log,vmrs_in,npro1,npro_max,
     >	press1_log,vmrs_out)
        do j=1,npro1
          vmrs1(j,i)=vmrs_out(j)
	  enddo
	enddo

c  ** get output filename **      
      print*,'Enter output file name (.prf or .ref):'
      read*,prffile_out

c  ** write output prf file **
	call zwriteprf(prffile_out,nvmr_max,npro_max,
     > iplanet,latitude,nvmr,npro1,molwt,idiso,height1,press1,temp1,
     > vmrs1)

	print*,'NB. Heights calculated by interpolation and'
      print*,'are probably not very accurate.'
      print*,'-Recalc using Profile'
      
      end


c-----------------------------------------------------------------------
	subroutine zreadprf(prffile,nvmr_max,npro_max,
     $   iplanet,latitude,nvmr,npro,molwt,idiso,height,press,temp,vmrs)
c-----------------------------------------------------------------------
c
c	read in a standard Nemesis .prf file
c
c	 variables
c	input:
c	 prffile		char	prf file
c	 nvmr_max		int	max no. of vmr profiles
c	 npro_max		int	max no. of pressure levels
c
c-----------------------------------------------------------------------
c	N. Teanby	10-02-04	Original code
c-----------------------------------------------------------------------

      implicit none
      integer nvmr_max,npro_max,iplanet,nvmr,npro
      integer iform,idiso(nvmr_max,2)
      real latitude,molwt,vmrs(npro_max,nvmr_max)
      real height(npro_max),press(npro_max),temp(npro_max)
      character*100 prffile
      character*1 buffer
      integer i,j,l
      
	open(10,file=prffile,status='old')
    
c  ** read comments and ignore **
1	continue
	read(10,*) buffer
      if (buffer.eq.'#') goto 1
      backspace(10)

c  ** read header **      
	read(10,*) iform
	read(10,*) iplanet,latitude,npro,nvmr,molwt
      do i=1,nvmr
        read(10,*) idiso(i,1),idiso(i,2)
      enddo
c  ** check max dimensions **
	if (nvmr.gt.nvmr_max) then
         print*, 'ERROR: zreadprf nvmr > nvmr_max'
         stop
        endif
	if (npro.gt.npro_max) then
         print*, 'ERROR: zreadprf npro > npro_max'
         stop
        endif
c  ** read vmrs **      
      read(10,*) buffer
      do i=1,npro
       read(10,*) height(i),press(i),temp(i),(vmrs(i,j),j=1,min(3,nvmr))
      enddo
	if (nvmr.gt.3) then
        do l=1,floor(real(nvmr-3-0.01)/6.)+1
          read(10,*) buffer
          do i=1,npro
            read(10,*) ( vmrs(i,j),j=4+6*(l-1),min(nvmr,3+(l*6)) )
          enddo
	  enddo
      endif
      close(10)
      
      return
      end

c-----------------------------------------------------------------------
	subroutine zlinint(n,np,x,y,nint,npint,xint,yint)
c-----------------------------------------------------------------------
c
c	linear interpolate n data points (x,y) to nint new x values xint
c
c	variables
c	input:
c	 n				int	number of points
c	 np				int	array dimensions
c	 x(np)			real	x values of series to interpolate
c	 y(np)			real	y values of series to interpolate
c	 nint				int	number of points in interpolated series
c	 npint			int	array dimensions of interpolated series
c	 xint(npint)		real	new x values to interpolate at
c	 y(np)			real	series to interpolate
c
c	output:
c	 yint(np)			real	interpolated series
c
c-----------------------------------------------------------------------
c	N. Teanby	15-10-03	Original code
c-----------------------------------------------------------------------

	implicit none
	integer i,j,n,np,nint,npint
	real y(np),yint(npint),x(np),xint(npint)


c  ** interpolation **
c  ** ascending order **
	if (x(1).lt.x(n)) then
       do i=1,nint
         if (xint(i).lt.x(1)) then
           yint(i)=y(1)
         else if (xint(i).gt.x(n)) then
           yint(i)=y(n)
         else
           do j=1,n-1
             if ((x(j).le.xint(i)).and.(x(j+1).ge.xint(i))) then
               yint(i) = y(j)+(y(j+1)-y(j))*(xint(i)-x(j))/(x(j+1)-x(j))
               goto 1
             endif
           enddo
           print*,'ERROR zlinint: interpolation error'
	   stop
1	     continue
         endif
       enddo
c  ** decending order **
	else
       do i=1,nint
         if (xint(i).gt.x(1)) then
           yint(i)=y(1)
         else if (xint(i).lt.x(n)) then
           yint(i)=y(n)
         else
           do j=1,n-1
             if ((x(j).ge.xint(i)).and.(x(j+1).le.xint(i))) then
               yint(i) = y(j)+(y(j+1)-y(j))*(xint(i)-x(j))/(x(j+1)-x(j))
               goto 2
             endif
           enddo
           print*, 'ERROR zlinint: interpolation error'
	   stop
2	     continue
         endif
       enddo
	endif
      
	return
	end
	
c-----------------------------------------------------------------------
	subroutine zwriteprf(prffile,nvmr_max,npro_max,
     $   iplanet,latitude,nvmr,npro,molwt,idiso,height,press,temp,vmrs)
c-----------------------------------------------------------------------
c
c	writen a standard Nemesis .prf file
c
c	 variables
c	input:
c	 prffile		char	prf file
c	 nvmr_max		int	max no. of vmr profiles
c	 npro_max		int	max no. of pressure levels
c
c-----------------------------------------------------------------------
c	N. Teanby	19-04-04	Original code
c-----------------------------------------------------------------------

      implicit none
      integer nvmr_max,npro_max,iplanet,nvmr,npro
      integer iform,idiso(nvmr_max,2)
      real latitude,molwt,vmrs(npro_max,nvmr_max)
      real height(npro_max),press(npro_max),temp(npro_max)
      character*100 prffile
      integer i,j,l
      
	open(20,file=prffile,status='unknown')
    
c  ** write header **      
	write(20,*) '           0'
	write(20,'(i3,2x,f8.2,2x,i5,2x,i5,2x,f8.3)')
     >	iplanet,latitude,npro,nvmr,molwt
      do i=1,nvmr
        write(20,'(i5,2x,i5)') idiso(i,1),idiso(i,2)
      enddo

c  ** write vmrs **      
      write(20,*) '  height (km)  press (atm)   temp (K)    vmrs....'
      do i=1,npro
        write(20,'(f12.3,2x,e12.5,2x,f12.4,3(2x,e12.5))')
     > height(i),press(i),temp(i),(vmrs(i,j),j=1,min(3,nvmr))
      enddo
	if (nvmr.gt.3) then
        do l=1,floor(real(nvmr-3-0.01)/6.)+1
          write(20,*) 'vmrs....'
          do i=1,npro
            write(20,'(e12.5,5(2x,e12.5))')
     >	( vmrs(i,j),j=4+6*(l-1),min(nvmr,3+(l*6)) )
          enddo
	  enddo
      endif
      close(20)
      
      return
      end
