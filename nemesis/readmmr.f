	subroutine readmmr(runname,l,n,k)
C	****************************************************************
C	Routine to read in complex refractive index data from a lookup
C	table at a single wavelength for a Maltmieser particle.	
C  
C	Input variable
C	runname	character*100	Name of file to read
C	l	real	Wavelength at which to search for complex ref index
C
C	Output variables
C 	n	real	Real refractive index value
C	k	real	Imaginary refractive index value
C
C	Ashwin Braude	14.12.2016	Original
C 
C     ****************************************************************

	implicit none

	character*100 runname, buffer, mm
	real n,k,l,lnew!,first,second,dlambda
	integer io, nspec, ispec
	logical mmexist
        integer idiag,iquiet
        common/diagnostic/idiag,iquiet

C	counthead = 1

	call file(runname,mm,'mmr')
	inquire(file=mm,exist=mmexist)
	if(mmexist)then
		open(46,file=mm,status='old')
56		read(46,1)buffer
1        	format(a)
         	if(buffer(1:1).eq.'#') goto 56!skip header
		read(46,*)nspec!read no of wavelengths in file
		read(46,*)ispec!if =1, data in wavelengths, if =0 data in wavenumbers
57		read(46,*,iostat=io)lnew,n,k
C2        	format(f4.2,f4.2,ES4.2)
		if(io.ne.0) then
			if(idiag.gt.0)print*,io
			goto 401!if eof reached, error message
		endif
		if(lnew.ne.l) goto 57
		close(46)
		goto 402!skip error message
	else
		print*,'Error readmmr: No reference table specified'
		print*,'for fixed complex refractive index spectrum' 
		print*,'in Coated Sphere model'
		stop
	endif

C	Error message
401	print*,'Error readmmr:'
	print*,'Input wavelength outside range of .mmr'
	stop

402	return

	end     
