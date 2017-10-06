	subroutine readrdw(runname,nconvfull,vconvfull,ngeom,nconv,
     1   vconv,rdwindices)
C	****************************************************************
C       Routine to read in a reduced wavelength grid
C  
C	Input variable
C	runname	character*100	Name of file to read
C 	nconvfull		integer		Number of wavelengths in full grid
C	vconvfull		real		Full wavelength grid array
C
C	Output variables
C 	nconv		integer		Number of wavelengths in reduced grid
C	vconv		real		Reduced wavelength grid array
C       rdwindices	integer		Subscripts of the full wavelength array where wavelengths in the reduced grid are found
C
C	Ashwin Braude	04.10.2017	Original
C 
C     ****************************************************************

	implicit none
        include '../radtran/includes/arrdef.f'
        INCLUDE 'arraylen.f'
	character*100 runname,buffer,rdw
	real tmp,vconv(mgeom,mconv),vconvfull(mgeom,mconv)
        integer rdwindices(mgeom,mconv)
	integer io,ispec,counter,cindex,cindex2,ngeom,igeom
        integer nconv(mgeom),nconvfull(mgeom),nconvtmp

        counter = 1!position in full wavelength grid
        cindex = 1!position in retrieved reduced grid
        cindex2 = 1!position in .rdw file

        print*,'Reading reduced wavelength grid .rdw'

	call file(runname,rdw,'rdw')
        do igeom=1,ngeom
	 open(46,file=rdw,status='old')
56	  read(46,1)buffer
1          format(a)
           if(buffer(1:1).eq.'#') goto 56!skip header

	  read(46,*)nconvtmp!read no of wavelengths in file
	  read(46,*)ispec!if =1, data in wavelengths, if =0 data in wavenumbers

          do while (cindex2.le.nconvtmp)
	   read(46,*,iostat=io)tmp
           if(counter.gt.nconvfull(igeom))then!if last wavelengths in reduced grid not present in full grid, skip last wavelengths
            goto 402
           endif
57         if(abs(tmp-vconvfull(igeom,counter)).lt.1e-5)then!if specified wavelength in reduced grid also present in full grid, store wavelength and subscript
            vconv(igeom,cindex)=tmp
            rdwindices(igeom,cindex)=counter
            counter=counter+1
            cindex=cindex+1
            cindex2=cindex2+1
           elseif(vconvfull(igeom,counter).gt.(1e-5+tmp))then!if specified wavelength in reduced grid not present in full grid, try again with next wavelength in reduced grid
            cindex2=cindex2+1
           else!if specified wavelength in full grid not found in reduced grid, try again with next wavelength in full grid
            counter=counter+1
            goto 57
           endif
          enddo

         close(46)

402      nconv(igeom)=cindex-1

         print*,'Reduced wavelength grid for igeom=',igeom
         print*,'has number of wavelengths ',nconv(igeom)
         print*,'Selected wavelengths:',vconv(igeom,1:nconv(igeom))

         counter=1!reset counters for next geometry
         cindex=1
         cindex2=1

        enddo

        return

	end     
