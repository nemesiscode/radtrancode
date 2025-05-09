      SUBROUTINE cirsrtfg_wave(runname,dist,inormal,iray,fwhm1,
     1 ispace,vwave,nwave,itype1, nem, vem, emissivity,tsurf,
     2 gradtsurf,nv,xmap,vconv,nconv,npath1,calcout,gradients,iscat)
C***********************************************************************
C_TITL:	CIRSRTFG_WAVE.f
C
C_DESC: 
C
C_ARGS:	Input variables:
C	runname		CHARA*100 Operation filename.
C	dist		REAL	Distance from Sun in units of AU.
C	inormal		INT	Flag for ortho:para ratio (0=equilibrium
C				=>1:1) (1=normal =>3:1).
C	fwhm1		REAL	Full-Width-at-Half-Max.
C       ispace          integer Indicates if wavelengths in vconv and
C				vwave are in wavenumbers(0) or 
C				wavelengths (1)
C	vwave(nwave)	REAL	Calculation wavenumbers.
C	nwave		INT	Number of calculation wavenumbers.
C	itype1		INT	Value designating the chosen scattering
C				routine (currently only scloud8 through
C				scloud11).
C	nv		INT	Number of 'on'/variable elements.
C	xmap(maxv,maxgas+2+maxcon,maxpro)	...
C			REAL	Matrix giving rate of change of elements
C				of T/P and aerosol .prf values with each
C				of the NV variables.
C	vconv(nconv)	REAL	Convolution wavenumbers.
C	nconv		INT	Number of convolution wavenumbers.
C	npath1		INT	Number of paths in calculation.
C
C       iscat           INT     Scattering flag mode
C
C	Output variables:
C	calcout(maxout3) REAL	Output values at each wavenumber for each
C				output type for each path
C	gradients(maxout4) REAL	Calculate rate of change of output with 
C				each of the NV variable elements, for
C				each wavenumber and for each path.
C
C_FILE:	unit=4		dump.out
C
C_CALL:	SUBPATHG	Reads in the .pat file, computes the atmospheric
C			absorber paths, then outputs the .drv file.
C	FILE		Forces file extension.
C	READ_KLIST	Reads in correlated-k files from .kls, passes
C			variables to common block INTERPK.
C	GET_SCATTER	Reads in scattering files (e.g. .sca).
C	GET_XSEC	Reads in cross-sections (.xsc) file.
C	CIRSRADG_WAVE
C	MAP2PRO		Converts the rate of change of output with respect
C			to layer properties to the rate of change of
C			output with .prf properties.
C	MAP2XVEC	Converts from rate of change of radiance with
C			profile .prf properties to user defined variables.
C	CIRSCONV	Convolves input spectrum with a bin of width
C			fwhm to produce an output spectrum.
C	CLOSE_SCAT	Closes any opened scattering files before next
C			iteration of retrieval algorithm.
C
C_HIST: 
C	30.7.2001 PGJI  Serious modification to deal with gradients.
C	7aug03	NT	corrected: inormal, real->integer and 
C			fwhm1, integer -> real. so consistent with other routines
C	29.2.2012 PGJI	Updated for Radtrans2.0
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
C     ../includes/arrdef.f defines the maximum values for a series of variables
C     (layers, bins, paths, etc.)
      INCLUDE '../includes/pathcom.f'
C     ../includes/pathcom.f holds the variables used by the software when 
C     calculating atmospheric paths (e.g. FWHM, IMOD, NPATH and LINKEY).
      INCLUDE '../includes/laycom.f'
C     ../includes/laycom.f holds variables used only by the path software
C     parameters are passed between routines mostly using common blocks
C     because of the extensive use of large arrays. NOTE: laycom uses
C     parameters defined in pathcom.
      INCLUDE '../includes/laygrad.f'
C     ../includes/laygrad.f holds the variables for use in gradient 
C     calculations.

      integer :: iparam,ipro,ntab1,ntab2
C     NTAB1: = NPATH*NWAVE, must be less than maxout3
C     NTAB2: = NPATH*NWAVE*NV, must be less than maxout4.
      integer :: i,j,nparam,jj,icont,status

      integer :: ilbl,ishape
      real :: gtsurf(maxpat)
      real :: RADIUS1
C     NB: The variables above have the added '1' to differentiate the 
C     variables passed into this code from that defined in
C     ../includes/pathcom.f. The definitions are explained above.

      real :: xref,xcomp,vv
C     XREF: % complete
C     XCOMP: % complete printed in increments of 10.
      character (len=100) :: drvfil,radfile,xscfil,runname,sfile
      character (len=100) :: klist, FWHMFILE
      integer :: iwave,ipath,k,igas,ioff1,ioff2,iv,ii,nkl
      real :: zheight(maxpro),radextra
      integer :: nsw,iswitch,rdamform,npath1
      logical :: scatterf,dustf,solexist,fexist,qfla,fwhmexist
      integer :: NFWHM,MFWHM
      parameter (MFWHM=1000)
      real :: VFWHM(MFWHM),XFWHM(MFWHM)

c     Allocatable arrays
      real, allocatable :: output(:),doutputdq(:,:,:)
      real, allocatable :: doutmoddq(:,:,:),doutdx(:,:),tempgtsurf(:)
      real, allocatable :: tempout(:),tgrad(:)
      real, allocatable :: y(:),yout(:),ygt(:),youtgt(:)    
      integer, allocatable :: isw(:)

c      real :: output(maxpat),doutputdq(maxpat,maxlay,maxgas+2+maxcon)
c      real ::  doutmoddq(maxpat,maxgas+2+maxcon,maxpro)
c      real :: doutdx(maxpat,maxv),tempgtsurf(maxout3)
c      real :: tempout(maxout3),tgrad(maxout4)
c      real :: y(maxout),yout(maxout),ygt(maxout),youtgt(maxout) 
c      integer :: isw(maxgas+2+maxcon)

C     Inputs and outputs
      real, intent(in) :: dist,fwhm1,tsurf
      integer, intent(in) :: inormal,iray,ispace,nwave,nconv
      integer, intent(in) :: itype1,nem,nv,iscat
      real, intent(in) :: vwave(nwave),vconv(nconv),vem(maxsec)
      real, intent(in) :: xmap(maxv,maxgas+2+maxcon,maxpro)
      real, intent(In) :: emissivity(maxsec)

      real, intent(out) :: gradtsurf(maxout3)
      real, intent(out) :: calcout(maxout3),gradients(maxout4)


C     Need simple way of passing planetary radius to nemesis
      INCLUDE '../includes/planrad.f'
      common/lbltable/ilbl

      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

c  ** non-LTE things - passed to cirsradg_wave using a common block **
	character*100 nltefile
	logical nlteexist
	integer inlte_flag
	real nlte_k1,nlte_A
      common/nlte_flags/inlte_flag,nlte_k1,nlte_A


C     ************************* CODE ***********************

C      PRINT*,'CIRSRTFG_WAVE - DEBUG'
C      PRINT*,'runname : ',runname
C      print*,'dist,inormal,iray,fwhm1,ispace',dist,inormal,iray,
C     1 fwhm1,ispace
C      print*,'nwave = ',nwave
C      do i=1,nwave
C       print*,i,vwave(i)
C      enddo
C      print*,'itype1,nem = ',itype1,nem
C      do i=1,nem
C       print*,i,vem(i),emissivity(i)
C      enddo
C      print*,'tsurf,nv,npath1,nconv = ',tsurf,nv,npath1,nconv
C      do i=1,nconv
C       print*,vconv(i)
C      enddo    


      if(idiag.gt.0)print*,'CIRSRTFG_WAVE - ILBL = ',ILBL
      IF(ILBL.EQ.2)THEN
         call file(runname,sfile,'sha')
         open(13,file=sfile,status='old')
         READ(13,*)ISHAPE
         close(13)

         if(idiag.gt.0)print*,'ISHAPE = ',ISHAPE
      ENDIF

C     Copy required FWHM to common block variable
      FWHM = FWHM1

C     See if file is present forcing FWHM to vary with wavelength/wavenumber
      CALL FILE(runname,FWHMFILE,'fwh')
      INQUIRE(FILE=FWHMFILE,EXIST=FWHMEXIST)
C     If such a file exists then read in the data
      IF(FWHMEXIST)THEN
         if(idiag.gt.0)then
          print*,'Reading FWHM information from : ',FWHMFILE
         endif
         OPEN(13,FILE=FWHMFILE,status='old')
          READ(13,*)NFWHM
          DO I=1,NFWHM
           READ(13,*)VFWHM(I),XFWHM(I)
          ENDDO
         CLOSE(13)
      ENDIF

c  ** set if an non-LTE .lte definition file is present **
      call file(runname,nltefile,'lte')
      inquire(file=nltefile,exist=nlteexist)
      if (nlteexist) then
         open(14,file=nltefile,status='old')
         read(14,*) inlte_flag
         if (inlte_flag.eq.1) then
c        * option 1 need k1 and A parameters *    
           read(14,*) nlte_k1,nlte_A
         endif
         close(14)
      else
         inlte_flag = 0
      endif

C Call subpathg to create layers, paths and the driver file.
C	  MOLWTX = XXMOLWT
      CALL subpathg(runname)
      npath1 = npath           ! npath is initilised in subpathg, set to
                               ! npath1 here so that it can be passed out

C Read the ktables or lbltables

      IF(ILBL.EQ.0)THEN
         CALL file (runname, klist, 'kls')
         if(idiag.gt.0)then
          WRITE(*,1050)klist

          WRITE(*,*)'     CALLING read_klist'
         endif
         CALL read_klist (klist, ngas, idgas, isogas, nwave, vwave, nkl)
         if(idiag.gt.0)then
          WRITE(*,*)'     read_klist COMPLETE'
          WRITE(*,*)' '
         endif
      ELSE
         CALL file (runname, klist, 'lls')
         if(idiag.gt.0)then
          WRITE(*,1050)klist

          WRITE(*,*)'     CALLING read_klbllist'
         endif
         CALL read_klbllist (klist, ngas, idgas, isogas, nwave, vwave)
         if(idiag.gt.0)then
          WRITE(*,*)'     read_klist COMPLETE'
          WRITE(*,*)' '
         endif
      ENDIF


C Now read the scattering files if required.
      scatterf = .FALSE.
      DO I=1,npath
        IF(imod(I).EQ.15.OR.imod(I).EQ.16)scatterf = .TRUE.      
      ENDDO

      QFLA = .FALSE.
      IF(ISCAT.EQ.5)THEN
        QFLA = .TRUE.
      ENDIF

      IF(qfla)scatterf = .TRUE.	

      IF(scatterf)THEN
        CALL file(runname, radfile, 'sca')
        CALL get_scatter(radfile,ncont)
      ENDIF

C ... and the xsc files likewise.
      dustf = .FALSE.
      IF(ncont.GT.0)dustf = .TRUE.

      IF(dustf)THEN
        CALL file(runname, xscfil, 'xsc')
        IF(qfla)THEN
		CALL get_xsec_red(xscfil,ncont,vwave)
	ELSE
		CALL get_xsec(xscfil, ncont)
	ENDIF
      ENDIF

      if(idiag.gt.0)print*,'NPATH = ',NPATH

C=======================================================================
C
C	Call CIRSradg_wave.	
C
C=======================================================================

      allocate(isw(maxgas+2+maxcon),STAT=status)
      if(status.gt.0)then
       print*,'error in CIRSrtfg_wave.f'
       print*,'isw was not properly allocated'
       stop
      endif

      nparam = ngas + 1 + ncont
      IF(flagh2p.EQ.1)nparam = nparam + 1

C Assess for which parameters a dr/dx value is actually needed
      nsw = 0
      DO iparam=1,nparam
        iswitch = 0
        DO iv=1,nv
          DO ipro=1,npro
            IF(xmap(iv,iparam,ipro).NE.0.0)iswitch = 1
          ENDDO
        ENDDO
        IF(iswitch.EQ.1)THEN
          nsw = nsw + 1
          isw(nsw) = iparam
        ENDIF
      ENDDO

      ntab1 = nwave*npath
      ntab2 = nwave*npath*nv
      IF(ntab1.GT.maxout3)THEN
        WRITE(*,*)'CIRSRTFG_WAVE.f :: Error: nwave*npath > maxout3'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'nwave : ',nwave
        WRITE(*,*)'npath : ',npath
        WRITE(*,*)'nwave*npath = ',ntab1,' maxout = ',maxout3
        STOP
      ENDIF
      IF(ntab2.GT.maxout4)THEN
        WRITE(*,*)'CIRSRTFG_WAVE.f :: Error: nwave*npath*nv > maxout4'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'nwave : ',nwave
        WRITE(*,*)'npath : ',npath
        WRITE(*,*)'nv : ',nv
        WRITE(*,*)'nwave*npath*nv = ',ntab2,' maxout4 = ',maxout4
        STOP
      ENDIF

C Determine the % complete ...

      if(idiag.gt.0)then
       WRITE(*,*)' '
       WRITE(*,*)' Waiting for CIRSradg_wave etal to complete ...'
      endif
      xcomp = 0.0


C Allocating variables
      allocate(output(maxpat),tempgtsurf(maxout3),STAT=status)
      if(status.gt.0)then
       print*,'error in CIRSrtfg_wave.f'
       print*,'output and tempgtsurf were not properly allocated'
       stop
      endif

      allocate(tempout(maxout3),tgrad(maxout4),STAT=status)
      if(status.gt.0)then
       print*,'error in CIRSrtfg_wave.f'
       print*,'tempout and tgrad were not properly allocated'
       stop
      endif

      DO 1000 iwave=1,nwave
        if(nwave.gt.1)then
          xref = 100.0*FLOAT(iwave - 1)/FLOAT(nwave - 1)
        else
          xref = 0.0
        endif
        IF(xref.GE.xcomp)THEN
          if(idiag.gt.0)WRITE(*,*)' Percent Complete = ',xcomp
          xcomp = xcomp + 10.0
        ENDIF

C       See what format the .prf file is in case we have to account for 
C       the fact that at each level the sum of vmrs is 1.
        AMFORM=rdamform(runname)

        vv = vwave(iwave)

C       Pass radius of planet to cirsradg_wave
C       radius2 is radius held in planrad common block. Pass this to
C       cirsradg_wave, in case it's been updated.
        radius1=radius2

        radextra=0.
        if(ipzen.eq.2)then
C          need the current height profile
           call readprfheight(runname,npro,zheight)
           radextra=zheight(npro)
        endif
  
        allocate(doutputdq(maxpat,maxlay,maxgas+2+maxcon),STAT=status)
        if(status.gt.0)then
         print*,'error in CIRSrtfg_wave.f'
         print*,'doutputdq was not properly allocated'
         stop
        endif

	CALL cirsradg_wave (dist,inormal,iray,delh,nlayer,npath,ngas,
     1  press,temp,pp,amount,iwave,ispace,AMFORM,vv,nlayin,
     2  layinc,cont,scale,imod,idgas,isogas,emtemp,itype1,
     3  nem, vem, emissivity, tsurf, gtsurf, RADIUS1,
     4  flagh2p,hfp,radextra,nsw,isw,output,doutputdq)

C Convert from rates of change with respect to layer variables to rates
C of change of .prf profile variables.
        allocate(doutmoddq(maxpat,maxgas+2+maxcon,maxpro),STAT=status)
        if(status.gt.0)then
         print*,'error in CIRSrtfg_wave.f'
         print*,'doutmoddq was not properly allocated'
         stop
        endif
        CALL map2pro(nparam,doutputdq,doutmoddq)
        deallocate(doutputdq)

C Convert from rates of change of .prf profile variables to rates of
C change of desired variables.
        allocate(doutdx(maxpat,maxv),STAT=status)
        if(status.gt.0)then
         print*,'error in CIRSrtfg_wave.f'
         print*,'doutdx was not properly allocated'
         stop
        endif
        CALL map2xvec(nparam,npro,npath,nv,xmap,doutmoddq,doutdx)
        deallocate(doutmoddq)

        DO ipath=1,npath
          ioff1 = nwave*(ipath - 1) + iwave
          if(ioff1.gt.maxout3)then
            print*,'cirsrtfg: ioff1 > maxout3'
            stop
          endif
          tempout(ioff1) = output(ipath)

C         Store ROC of surface temperature too
          tempgtsurf(ioff1)=gtsurf(ipath)
C          print*,'YYY',iwave,ipath,ioff1,tempout(ioff1),
C     1      tempgtsurf(ioff1)

          DO iv=1,nv
            ioff2 = nwave*nv*(ipath - 1) + (iv - 1)*nwave + iwave
            if(ioff2.gt.maxout4)then
             print*,'cirsrtfg_wave: ioff2 > maxout4'
             stop
            endif
            tgrad(ioff2) = doutdx(ipath,iv)
          ENDDO
        ENDDO
        deallocate(doutdx)
1000  CONTINUE

      deallocate(output,isw)

C Convolve output calculation spectra with a square bin of width FWHM to
C get convoluted spectra.
      allocate(yout(maxout),youtgt(maxout),STAT=status)
      if(status.gt.0)then
       print*,'error in CIRSrtfg_wave.f'
       print*,'yout and youtgt were not properly allocated'
       stop
      endif

      allocate(y(maxout),ygt(maxout),STAT=status)
      if(status.gt.0)then
       print*,'error in CIRSrtfg_wave.f'
       print*,'y and ygt were not properly allocated'
       stop
      endif

      DO ipath=1,npath

        DO iwave=1,nwave
           ioff1 = nwave*(ipath - 1) + iwave
           y(iwave) = tempout(ioff1)
           ygt(iwave) = tempgtsurf(ioff1)
        ENDDO		


        if(Ilbl.eq.0)then
         CALL cirsconv(runname,fwhm,nwave,vwave,y,nconv,vconv,yout,
     1  FWHMEXIST,NFWHM,VFWHM,XFWHM)	
         CALL cirsconv(runname,fwhm,nwave,vwave,ygt,nconv,vconv,
     1	youtgt,FWHMEXIST,NFWHM,VFWHM,XFWHM)
        else
         CALL lblconv1(runname,fwhm,ishape,nwave,vwave,y,nconv,
     1     vconv,yout,FWHMEXIST,NFWHM,VFWHM,XFWHM)	
         CALL lblconv1(runname,fwhm,ishape,nwave,vwave,ygt,nconv,
     1     vconv,youtgt,FWHMEXIST,NFWHM,VFWHM,XFWHM)
        endif
        DO iconv=1,nconv
          ioff1 = nconv*(ipath - 1) + iconv
          calcout(ioff1) = yout(iconv)
          gradtsurf(ioff1) = youtgt(iconv)
C          print*,'ZZZ',iconv,ipath,ioff1,calcout(ioff1),
C     1     gradtsurf(ioff1)
        ENDDO

        DO iv=1,nv
          DO iwave=1,nwave
            ioff2 = nwave*nv*(ipath - 1) + (iv - 1)*nwave + iwave
            y(iwave) = tgrad(ioff2)
          ENDDO

          if(ilbl.eq.0)then
           CALL cirsconv(runname,fwhm,nwave,vwave,y,nconv,vconv,yout,
     1 FWHMEXIST,NFWHM,VFWHM,XFWHM)
          else
           CALL lblconv1(runname,fwhm,ishape,nwave,vwave,y,nconv,
     1      vconv,yout,FWHMEXIST,NFWHM,VFWHM,XFWHM)
          endif

          DO iconv=1,nconv
            ioff2 = nconv*nv*(ipath - 1) + (iv - 1)*nconv + iconv
            gradients(ioff2) = yout(iconv)
C            print*,'YYY',iconv,iv,ipath,ioff2,gradients(ioff2)
          ENDDO
        ENDDO
      ENDDO

      deallocate(tempgtsurf,tempout)
      deallocate(yout,youtgt)
      deallocate(y,ygt)

C=======================================================================
C
C	Wrap up: formats, return, and end.
C
C=======================================================================

      DO ii=101,nkl+100
        CLOSE(UNIT=ii)
      ENDDO

      CALL close_scat(ncont)

1050  FORMAT (/'Klist filename: ', A)

      RETURN

      END
