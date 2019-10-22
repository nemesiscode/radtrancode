	PROGRAM CIRSdrvg_wavePY
C*********************************************************************
C_TITL:	CIRSdrvg_wavePY
C
C_DESC:	Driver program for new gradient-version of cirsrad.
C
C_ARGS:	
C
C_CALL: 
C	file		Forces file extension.
C	remsp		Remove space from passed string.
C	upcase		Make upper-case the passed string.
C	cirsrtfg
C
C_HIST:	30/7/01	PGJI	ORIGINAL VERSION.
C       29/2/12 PGJI	Updated for Radtrans2.0

C*********************************************************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
C ../includes/arrdef.f defines the maximum values for a series of variables
C (layers, bins, paths, etc.)

      INTEGER intmod,nconv,nwave,i,j,itype,imie,imie1,iscat,ishape
      INTEGER iconv,ioff1,ioff2,iv,nx,iwave,ipath,ispace,lcdr
      INTEGER npoints,planet,inormal,iray,iptf,npath,ichannel,nem
      INTEGER k,ilbl

      REAL vconv(maxbin),vwave(maxbin),vkstart,vkend,vkstep
      REAL vref,delvk,vmin,vmax,vem(MAXSEC),emissivity(MAXSEC),tsurf
      REAL dist,wnumbot,wnumtop,delv,fwhm,gradtsurf(maxout3)
c      real :: calcout_noconv(maxout3),gradients_noconv(maxout4)

      REAL tot_time
      real rate
      integer c1,c2,cr,time1,time2,cm,npro,nvmr,ncont

      REAL calcout(maxout3),gradients(maxout4)
      REAL xmap(maxv,maxgas+2+maxcon,maxpro)
      REAL y(maxout),ygt(maxout4)
      CHARACTER*80 text
      CHARACTER*100 opfile,runname1,sfile,outfile,ciafil
      CHARACTER*1 ANS

      REAL xwave(100),xf(10,100),xg1(10,100),xg2(10,100)
      REAL tnco,twave,frac,tico
      COMMON /hgphas/ xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /imiescat/imie1
      COMMON /LBLTABLE/ ILBL

      INCLUDE '../includes/ciacom.f'
      INCLUDE '../includes/gascom.f'
      INCLUDE '../includes/planrad.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA

C     Solar spectrum variables and flags
      integer iform,iread,solnpt
      real solwave(maxbin),solrad(maxbin),solradius
      character*100 solfile,solname
      logical solexist
      common/solardat/iread, iform, solradius, solwave, solrad,  solnpt


C-----------------------------------------------------------------------
C
C	Write a header
C
C-----------------------------------------------------------------------
      jradf=-1
      jloggf=-1


C     Read in gas information
      CALL RESERVEGAS

C Obtain the system time for use in determining the total time used by the
C program since execution.
      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)

      CALL system_clock(time1)

C Reset xwave to force read of hgphase* files for scloud10-11.
      xwave(1)=-1

      WRITE(*,*)' '
      WRITE(*,*)'           WELCOME TO CIRSDRVG_WAVEPY'
      WRITE(*,*)' '

C-----------------------------------------------------------------------
C
C	Prompt the user for prescribed variables.
C
C-----------------------------------------------------------------------

      WRITE(*,*)'Give operation filename: '
      READ(*,1000)opfile
1     FORMAT(A)


C-----------------------------------------------------------------------
C
C       Reading all required parameters by CIRSrtf_wave from .cdr file
C
C-----------------------------------------------------------------------

      lcdr=37
      CALL file(opfile,runname1,'gdr')
      open(lcdr,file=runname1,status='old')
      print*,runname1
      read(lcdr,*)dist
      read(lcdr,*)fwhm
      read(lcdr,*)ispace
      read(lcdr,*)ilbl
      read(lcdr,*)itype
      read(lcdr,*)npath
      read(lcdr,*)nconv
      do i=1,nconv
        read(lcdr,*)vconv(i)
      enddo
      read(lcdr,*)nem
      do i=1,nem
        read(lcdr,*)vem(i),emissivity(i)
      enddo
      read(lcdr,*)tsurf
      read(lcdr,*)nx
      read(lcdr,*)nvmr
      read(lcdr,*)ncont
      read(lcdr,*)npro
      do i=1,nx
       do j=1,nvmr+3+ncont
        do k=1,npro 
         read(lcdr,*)xmap(i,j,k)
        enddo
       enddo
      enddo
      close(lcdr)

      print*,'nconv = ',nconv


C-----------------------------------------------------------------------
C
C       Calculating calculation wavenumbers/wavelengths
C
C-----------------------------------------------------------------------

      if(ilbl.eq.2)then
        CALL readkklblhead(opfile,vkstart,vkend,vkstep)
        if(fwhm.gt.0.0)then  !if it is lower then we read .fil file
         call file(opfile,sfile,'sha')
         open(13,file=sfile,status='old')
         READ(13,*)ISHAPE
         close(13)
        endif
        print*,'fwhm = ',fwhm
        print*,'ishape = ',ishape
        call wavesetc(opfile,vkstart,vkend,vkstep,nconv,vconv,
     1   fwhm,ishape,nwave,vwave)
        print*,'nwave = ',nwave
c        do i=1,nwave
c         print*,vwave(i)
c        enddo 
c        pause
      endif


      if(ilbl.eq.0)then
        CALL readkkhead(opfile,vkstart,vkend,vkstep)
        print*,'fwhm = ',fwhm
        call wavesetb(opfile,vkstart,vkend,vkstep,nconv,vconv,
     1   fwhm,nwave,vwave)
        print*,'nwave = ',nwave
      endif





C-----------------------------------------------------------------------
C
C       Get scattering type and Ortho/Para information.
C
C-----------------------------------------------------------------------

      CALL READFLAGS(OPFILE,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1 IPTF,IMIE,IUVSCAT)
      IMIE1=IMIE

      CALL FILE(opfile,CIAFIL,'cia')

      OPEN(12,FILE=CIAFIL,STATUS='OLD')
       READ(12,1)ANAME
       READ(12,*) DNU
       READ(12,*) IPARA
      CLOSE(12)
      IREAD1=1
      IREAD2=1
      IF(IPARA.EQ.0)THEN
       ANAME1=ANAME   
       DNU1=DNU
      ELSE
       ANAME2=ANAME     
       DNU2=DNU  
       IPARA2=IPARA
      ENDIF


C-----------------------------------------------------------------------
C
C	Call CIRSrtfg_wave
C
C-----------------------------------------------------------------------
	


C     See if there is a solar or stellar reference spectrum and read in
C     if present.
      call file(opfile,solfile,'sol')
      inquire(file=solfile,exist=solexist)
      if(solexist)then
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
      endif
      iform=0

      CALL cirsrtfg_wave (opfile, dist, inormal, iray, fwhm, ispace, 
     1 vwave,nwave,itype, nem, vem, emissivity, tsurf, gradtsurf, nx, 
     1 xmap, vconv, nconv, npath, calcout, gradients,iscat)

      CALL FILE(opfile,outfile,'out')
      OPEN(2,FILE=outfile,STATUS='unknown')
      WRITE(2,*)npath,nconv,nx
      DO ipath=1,npath
	WRITE(2,*)ipath
        DO iconv=1,nconv
          ioff1 = nconv*(ipath-1) + iconv
          WRITE(2,*)vconv(iconv),calcout(ioff1)
        ENDDO
        DO iv=1,nx
          DO iconv=1,nconv
            ioff2 = nconv*nx*(ipath - 1) + (iv - 1)*nconv + iconv
            y(iconv)=gradients(ioff2)
          ENDDO
          DO iconv=1,nconv
           WRITE(2,*)y(iconv)
          ENDDO
c          WRITE(2,*)(y(iconv),iconv=1,nconv)
        ENDDO
      ENDDO
      CLOSE(2)

c      CALL FILE(opfile,outfile,'gut')
c      OPEN(3,FILE=outfile,STATUS='unknown')
c      WRITE(3,*)npath,nwave,nx
c      DO ipath=1,npath
c        WRITE(3,*)ipath
c        DO iwave=1,nwave
c          ioff1 = nwave*(ipath-1) + iwave
c          WRITE(3,*)vwave(iwave),calcout_noconv(ioff1)
c        ENDDO
c        DO iv=1,nx
c          DO iwave=1,nwave
c            ioff2 = nwave*nx*(ipath - 1) + (iv - 1)*nwave + iwave
c            ygt(iwave)=gradients_noconv(ioff2)
c          ENDDO
c          DO iwave=1,nwave
c           WRITE(3,*)ygt(iwave)
c          ENDDO
c        ENDDO
c      ENDDO
c      CLOSE(3)

C----------------------------------------------------------------------- 
C 
C	Wrap up: formats and end.
C
C-----------------------------------------------------------------------

      WRITE(*,*)' CIRSDRVG_WAVE.f :: calculation complete.'
      CALL system_clock(time2)
      tot_time =  (time2-time1)/rate
      WRITE(*,200)tot_time
200   FORMAT(' Elapsed time including convolution (sec)= ',F8.1)

1000  FORMAT (A60)

      END
