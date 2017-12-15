	subroutine readrdw(runname,nconvfull,vconvfull,ngeom,nconv,
     1   vconv,rdwindices,rank,irank,nconvi,vconvi,rdwindicesi,
     2   maxirank,rankdiff)
C	****************************************************************
C       Routine to read in a reduced wavelength grid, and rank each wavelength so 
C       that the most important wavelengths are used in the earliest iterations of
C	Nemesis, while lesser-ranked wavelengths are saved to later iterations.
C  
C	Input variable
C	runname	character*100	Name of file to read
C 	nconvfull		integer		Number of wavelengths in full grid
C	vconvfull		real		Full wavelength grid array
C	irank			irank		Highest wavelength rank (lowest in value) to be taken into account in given set of iterations
C
C	Output variables
C 	nconv		integer		Number of wavelengths in reduced grid
C	vconv		real		Reduced wavelength grid array
C       rdwindices	integer		Subscripts of the full wavelength array where wavelengths in the reduced grid are found
C       rank		integer		'Ranking' of each wavelength
C 	nconvi		integer		Number of wavelengths of specific rank in reduced grid
C	vconvi		real		Reduced wavelength grid array (specific rank only)
C       rdwindicesi	integer		Subscripts of the full wavelength array where wavelengths of specific rank in the reduced grid are found
C       maxirank        integer         Maximum rank displayed in the .rdw file
C       rankdiff        integer         Difference between old rank and new rank (usually =1 unless towards end of retrieval)
C
C	Ashwin Braude	04.10.2017	Original
C 
C     ****************************************************************

	implicit none
        include '../radtran/includes/arrdef.f'
        INCLUDE 'arraylen.f'
	character*100 runname,buffer,rdw
	real tmp1,tmp2,vconv(mgeom,mconv),vconvfull(mgeom,mconv)
        integer rdwindices(mgeom,mconv),rank(mgeom,mconv)
	integer io,ispec,counter,cindex,cindex2,ngeom,igeom,irank
        integer nconv(mgeom),nconvfull(mgeom),nconvtmp,nconvi(mgeom)
        integer cindex3,rdwindicesi(mgeom,mconv),maxirank,rankdiff
        real vconvi(mgeom,mconv)

        counter = 1!position in full wavelength grid
        cindex = 1!position in retrieved reduced grid
        cindex2 = 1!position in .rdw file
        cindex3 = 1!position in wavelength grid containing only wavelengths of irank

        print*,'Reading reduced wavelength grid .rdw'

	call file(runname,rdw,'rdw')
        do igeom=1,ngeom
	 open(46,file=rdw,status='old')
56	  read(46,1)buffer
1          format(a)
           if(buffer(1:1).eq.'#') goto 56!skip header

	  read(46,*)nconvtmp!read no of wavelengths in file
	  read(46,*)ispec!if =1, data in wavelengths, if =0 data in wavenumbers
          read(46,*)maxirank!maximum ranking of wavelengths displayed in .rdw file
          if(maxirank.lt.1)then
           print*,'Error readrdw: Maxirank must be greater than 0'
           stop
          endif

          do while (cindex2.le.nconvtmp)
	   read(46,*,iostat=io)tmp1,tmp2
           if(counter.gt.nconvfull(igeom))then!if last wavelengths in reduced grid not present in full grid, skip last wavelengths
            goto 402
           endif
57         if(abs(tmp1-vconvfull(igeom,counter)).lt.1e-5)then!if specified wavelength in reduced grid also present in full grid, store wavelength and subscript
            if(tmp2.le.irank)then!only select wavelengths with correct rank
             vconv(igeom,cindex)=tmp1
             rank(igeom,cindex)=tmp2
             rdwindices(igeom,cindex)=counter
             if((tmp2.le.irank).and.(tmp2.gt.irank-rankdiff))then
              vconvi(igeom,cindex3)=tmp1
              rdwindicesi(igeom,cindex3)=cindex
              cindex3 = cindex3 + 1
             endif
             counter=counter+1
             cindex=cindex+1
             cindex2=cindex2+1
            else!if a wavelength is of the wrong rank, go to next wavelength in reduced grid
             counter=counter+1
             cindex2=cindex2+1
            endif
           elseif(vconvfull(igeom,counter).gt.(1e-5+tmp1))then!if specified wavelength in reduced grid not present in full grid, try again with next wavelength in reduced grid
            cindex2=cindex2+1
           else!if specified wavelength in full grid not found in reduced grid, try again with next wavelength in full grid
            counter=counter+1
            goto 57
           endif
          enddo

         close(46)

402      nconv(igeom)=cindex-1
         nconvi(igeom)=cindex3-1

         if(irank.le.maxirank)then
          print*,'Reduced wavelength grid for igeom=',igeom
          print*,'has number of wavelengths ',nconv(igeom)
          print*,'Selected wavelengths:',vconv(igeom,1:nconv(igeom))
         endif

         counter=1!reset counters for next geometry
         cindex=1
         cindex2=1
         cindex3=1

        enddo

        return

	end     
