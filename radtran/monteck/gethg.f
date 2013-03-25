      subroutine gethg(opfile,ncont,nlambda,xlambda,phased)
C     $Id: gethg.f,v 1.1 2005/08/01 11:11:57 irwin Exp $
C     ****************************************************************
C
C     Subroutine to read in H-G scattering parameter spectra from a 
C     .xsc and hgphase*.dat files
C
C
C     Input variables
C	opfile		character*100	name of x-section file
C
C     Output variables
C       ncont	integer	Number of particle types
C	nlambda	integer	Number of wavenumbers/wavelengths for scattering 
C                       property spectra
C	xlambda(nlambda) real Wavenumbers/wavelengths
C 	phased(maxcon,maxsec,5) real scattering spectra:
C				     xsection,omega,f,g1,g2
C
C     Revised  Pat Irwin  1/8/05
C     
****************************************************************

      implicit none
      include '../includes/arrdef.f'
      character*100 opfile,phafil,xscfil,buffer
      integer nlambda,i,j,k,ncont,iunit
      real xlambda(maxsec),phased(maxcon,maxsec,5)
      real xsec(2,maxcon,maxsec),vsec(maxsec),omega(maxcon)
C     ****************************************************************



      call file(opfile,xscfil,'xsc')
      iunit=19
      open(iunit,file=xscfil,status='old')
1     format(a)
54    read(iunit,1)buffer
      if(buffer(1:1).eq.'#')goto 54
      read(buffer,*)ncont

      j=0
105   j=j+1
      if(j.gt.maxsec)then
       print*,'Error in get_hg: Too many wavelengths in xsc file'
       print*,'j, maxsec : ',j,maxsec
       stop
      endif
      read(iunit,*,end=106)vsec(j),(xsec(1,i,j),i=1,ncont)
      read(iunit,*,end=106)(omega(i),i=1,ncont)
      do 344 i=1,ncont
C       xsec(2,i,j)=omega(i)*xsec(1,i,j)
       xsec(2,i,j)=omega(i)
344   continue
      goto 105
106   continue
      nlambda=j-1

      close(iunit)
      print*,'Xsc file read OK'
      print*,'Ncont is', ncont
      print*,' '
      print*,'Now reading hgphase*.dat files...'
C     Now read in hgphase*.dat files
      do 90 j=1,ncont
        phafil='hgphase*.dat'
        phafil(8:8) = char(j+48)
        print*,phafil
        open(iunit,file=phafil,status='old')
        do 110 i=1,nlambda
         read(iunit,*)xlambda(i),(phased(j,i,k),k=3,5)
         phased(j,i,1)=xsec(1,j,i)
         phased(j,i,2)=xsec(2,j,i)
         print*,xlambda(i),(phased(j,i,k),k=1,5)        
110     continue
        close(iunit)
90     continue

      return

      end
