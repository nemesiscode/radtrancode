      subroutine dumpflux(LUNIS,ioff,nlays,nmu,nf,radg,umif,uplf,iwave,
     1  x,nwave,Ig,ng,j0,xmu,solar,taus)
C     ***************************************************************
C     Subroutine to format output of scloud11flux into a more usable
C     internal radiation field and output to a meta file.
C
C     Input variables
C	LUNIS	integer	Unit number of output metafile
C	ioff	integer	offset record number in output metafile
C	nlays	integer	Number of scattering angles
C	nmu	integer	Number of zenith quadrature angles
C	nf	integer	Number of fourier components
C	radg(maxmu)	real	Ground radiance field
C	umif(maxmu,maxscatlay,maxf) real umif Radiance field going
C                                        upwards out of the top of each
C					 layer
C	uplf(maxmu,maxscatlay,maxf) real umif Radiance field going
C                                        downwards out of the bottom of each
C					 layer
C	iwave	integer	Wavelength ordinate in metafile
C	x	real	Wavelength/wavenumber No longer used.
C	nwave	integer	Number of wavelengths
C	Ig	integer	Ordinate of k-distribution integration
C	ng	integer	Number if g-ordinates
C	j0	integer	Nearest quadrature zenith angle to solar
C			zenith angle.
C	xmu	real	Cos(zolar zenith angle)
C	solar	real	Solar flux at this wavelength (weighted 
C			by 2*!pi*wt1(j0)
C	taus(maxscatlayer)	real	Optical depth of each layer
C
C     Documented by Pat Irwin	13/3/14
C
C     ***************************************************************
      implicit none
      INCLUDE '../includes/arrdef.f'
      integer LUNIS,ioff,nwave,irec,j0
      integer nlays,nmu,iwave,Ig,ng,joff,i,j,k,i1,nf
      real x,umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real radg(maxmu),tmpzero,taus(maxscatlay),xmu,solar,tmp,sol1

      irec = ioff+2*((iwave-1)*ng + (ig-1))*nlays*nmu*(nf+1)

C     Add in bottom radiance to complete radiance field at base of
C     atmosphere. NB. Scattering flux defines uplf to be radiance going 
C     downwards out of bottom of layer and umif to be the radiance going
C     upwards out of the top. We change the meaning here to be the 
C     total radiance field (up and down) at the bottom of each layer.  
C
C     In addition, the order of layers of the umif,uplf arrays is from 
C     top to bottom. We here reverse the order to be from bottom to top
C     so that the alignment with the baseh array is correct.

C     Also deduct direct solar part here as it screws up the fourier
C     decomposition

C     First strip out direct solar component from uplf
      tmp=0.0
      do j=1,nlays
       tmp = tmp+taus(j)/xmu
       sol1 = solar*exp(-tmp)

       do i1=1,nf+1
        uplf(j0,j,i1)=uplf(j0,j,i1)-sol1
       enddo


      enddo

      tmpzero=0.0
C     Reverse order of layers here. Scattering code has 1 as the top
C     layer. Rest of code have 1 as the bottom layer.
 
      do k=1,nmu
        do i1=1,nf+1
          if(i1.eq.1) then
            write(LUNIS,REC=IREC)radg(k)
          else
            write(LUNIS,REC=IREC)tmpzero
          endif
          irec=irec+1
          write(LUNIS,REC=IREC)uplf(k,nlays,i1)
          irec=irec+1
        enddo
      enddo
      
      do j=nlays-1,1,-1
       do k=1,nmu
        do i1=1,nf+1
         write(LUNIS,REC=IREC)umif(k,j+1,i1)
         irec=irec+1
         write(LUNIS,REC=IREC)uplf(k,j,i1)
         irec=irec+1
        enddo
       enddo
      enddo



      return

      end
