      subroutine wavesetc(runname,vkstart,vkend,vkstep,nconv,vconv,
     1 fwhm,ishape,nwave,vwave)
C     $Id:
C     ******************************************************************
C     Subroutine to calculate which 'calculation' wavelengths are needed
C     to cover the required 'convolution wavelengths'.
C
C     Input variables
C	runname	character*100 Run-name
C	vkstart	real	start wavenumber of k-tables
C	vkend	real	end wavenumber of k-tables
C	vkstep  real	wavenumber step of ktables
C	nconv 	integer	Number of convolution wavelengths
C	vconv(mconv) real	Convolution wavenumbers
C	fwhm	real	FWHM of convolved spectrum
C
C     Output variables
C	nwave	integer Number of calculation wavenumbers
C	vwave(mwave) real Calculation wavenumbers
C 
C
C     Pat Irwin	19/5/97		Original
C     Pat Irwin 21/3/00		Converted from NIMS
C     Pat Irwin 17/10/03	Tidied for Nemesis
C
C     ******************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer nconv,i,j,k,ico,nco,nwave,j1,j2,nconv1,nsub,jj
      integer ishape
      real vwave(mwave),v1,v2,save(100000),temp,vconv(mconv),vcentral
      real xdiff,test,vkstart,vkend,vkstep,fwhm,vfil(1000),dv,vj
      real fil(1000)
      logical flag
      character*100 runname
C     *******************************************************************

C     If k-tables are tabulated on an irregular wavelength grid OR have
C     have been calculated with a built-in instrument function then
C     just set the calculation wavelengths to be the convolution wavelengths


C      print*,'vconv'
C      print*,nconv
C      do i=1,nconv
C       print*,i,vconv(i)
C      enddo

      ico = 0
      if(fwhm.lt.0.0)then
        call file(runname,runname,'fil')
        print*,'Reading : ',runname
        open(12,file=runname,status='old')
         read(12,*)nconv1
         if(nconv.ne.nconv1)then
          print*,'Error in wavesetc'
          print*,'Number of channels in filter file not the same as'
          print*,'number of channels is spectral measurement file'
          print*,'nconv1,nconv : ',nconv1,nconv
          stop
         endif
         do 442 k=1,nconv1
          read(12,*)vcentral
          read(12,*)nsub
          if(nsub.gt.1000)then
           print*,'Error in wavesetc: vfil array not big enough'
           print*,nsub
           stop
          endif
          do j=1,nsub
           read(12,*)vfil(j),fil(j)
          enddo

          do 443 i=1,nconv
           dv = abs(vcentral-vconv(i))

           if(dv.lt.0.00001)then
            j1=((vfil(1)-vkstart)/vkstep)-1
            j2=((vfil(nsub)-vkstart)/vkstep)+1
            v1 = vkstart + (j1-1)*vkstep
            v2 = vkstart + (j2-1)*vkstep
            if(v1.lt.vkstart.or.v2.gt.vkend)then
              print*,'Warning from wavesetc'
              print*,'Channel wavelengths not covered by ktables'
              print*,'v1,v2,vkstart,vkend',v1,v2,vkstart,vkend
            endif

            do jj=j1,j2
             ico=ico+1
             vj = vkstart + jj*vkstep
             save(ico)=vj
            enddo

           endif

443       continue

442      continue

        close(12)

      else

        do 444 i=1,nconv
          if(ishape.eq.0)then
           dv = 0.5*fwhm
          elseif(ishape.eq.1)then
           dv = fwhm
          elseif(ishape.eq.2)then
           dv = 3*0.5*fwhm/sqrt(alog(2.0))
          else
           dv = 3*fwhm
          endif
          j1=int((vconv(i)-dv-vkstart)/vkstep)
          j2=2+int((vconv(i)+dv-vkstart)/vkstep)
          v1 = vkstart + (j1-1)*vkstep
          v2 = vkstart + (j2-1)*vkstep
          if(v1.lt.vkstart.or.v2.gt.vkend)then
              print*,'Warning from wavesetc'
              print*,'Channel wavelengths not covered by ktables'
              print*,'v1,v2,vkstart,vkend',v1,v2,vkstart,vkend
          endif

          do j=j1,j2
           v2 = v1 + (j-j1)*vkstep
           if(v2.ge.vkstart.and.v2.le.vkend)then
            ico = ico+1
            save(ico)=v2
           endif
          end do

444     continue

      endif

      nco=ico

C     sort calculation wavelengths into order
440   continue
      flag=.false.
      do 445 i=1,nco-1
        if(save(i).gt.save(i+1))then
         flag=.true.
         temp = save(i)
         save(i)=save(i+1)
         save(i+1)=temp
        endif
445   continue
      if (flag) goto 440
      nwave=nco
      do i=1,nwave
       vwave(i)=save(i)
      enddo

C      print*,'nco = ',nco
C     Now weed out repeated wavelengths
      vwave(1)=save(1)
      ico=1
      xdiff = 0.9*vkstep
      if(xdiff.gt.0)then
      print*,'Weeding out repeated wavelengths. vkstep = ',
     &		vkstep
      do 446 i=2,nco
       test = abs(save(i)-vwave(ico))
       if(test.ge.xdiff)then
         ico=ico+1
         if(ico.gt.mwave)then
          print*,'Overflow in wavesetc.f'
          print*,'ico > mwave',ico,mwave
          print*,save(i)
          stop
         endif
         vwave(ico)=save(i)
C         print*,'ico,vwave',ico,vwave(ico)
       end if
446   continue
      nwave=ico
      endif

C      print*,'nwave = ',nwave
C      do i=1,nwave
C       print*,i,vwave(i)
C      enddo

      print*,'wavesetc: nwave = ',nwave
      return
      end
