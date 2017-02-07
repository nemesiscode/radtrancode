      SUBROUTINE cirsrtfg_wave(runname,dist,inormal,iray,fwhm1,
     1 ispace,vwave,nwave,itype1, nem, vem, emissivity,tsurf,
     2 gradtsurf,nv,xmap,vconv,nconv,npath1,calcout,gradients)
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

      INTEGER iparam,ipro,ntab1,ntab2
C     NTAB1: = NPATH*NWAVE, must be less than maxout3
C     NTAB2: = NPATH*NWAVE*NV, must be less than maxout4.
      INTEGER nwave,i,j,nparam,jj,icont,nconv

      INTEGER itype1,npath1,ispace,nem,ILBL,ishape
      REAL fwhm1,vem(maxsec),emissivity(maxsec),tsurf
      real gtsurf(maxpat)
      REAL RADIUS1
C     NB: The variables above have the added '1' to differentiate the 
C     variables passed into this code from that defined in
C     ../includes/pathcom.f. The definitions are explained above.

      REAL dist,xref,xcomp
      INTEGER inormal,iray
C     XREF: % complete
C     XCOMP: % complete printed in increments of 10.
      REAL vwave(nwave),output(maxpat),vconv(nconv)
      REAL doutputdq(maxpat,maxlay,maxgas+2+maxcon)
      REAL doutmoddq(maxpat,maxgas+2+maxcon,maxpro)
      REAL doutdx(maxpat,maxv),gradtsurf(maxout3)
      REAL tempgtsurf(maxout3)
      REAL tempout(maxout3),tgrad(maxout4)
      REAL calcout(maxout3),gradients(maxout4)
      REAL xmap(maxv,maxgas+2+maxcon,maxpro)
      REAL y(maxout),yout(maxout),vv,ygt(maxout),youtgt(maxout)
      CHARACTER*100 drvfil,radfile,xscfil,runname,sfile
      CHARACTER*100 klist
      INTEGER iwave,ipath,k,igas,ioff1,ioff2,iv,nv
      REAL zheight(maxpro),radextra
      INTEGER nsw,isw(maxgas+2+maxcon),iswitch,rdamform
      LOGICAL scatterf,dustf,solexist,fexist

C     Need simple way of passing planetary radius to nemesis
      INCLUDE '../includes/planrad.f'
      common/lbltable/ilbl

C     ************************* CODE ***********************

      PRINT*,'CIRSRTF_WAVE - ILBL = ',ILBL
      IF(ILBL.EQ.2)THEN
         call file(runname,sfile,'sha')
         open(13,file=sfile,status='old')
         READ(13,*)ISHAPE
         close(13)

         print*,'ISHAPE = ',ISHAPE
      ENDIF


C Call subpathg to create layers, paths and the driver file.
C	  MOLWTX = XXMOLWT
      CALL subpathg(runname)
      npath1 = npath           ! npath is initilised in subpathg, set to
                               ! npath1 here so that it can be passed out

C Read the ktables or lbltables

      IF(ILBL.EQ.0)THEN
         CALL file (opfile, klist, 'kls')
         WRITE(*,1050)klist

         WRITE(*,*)'     CALLING read_klist'
         CALL read_klist (klist, ngas, idgas, isogas, nwave, vwave)
         WRITE(*,*)'     read_klist COMPLETE'
         WRITE(*,*)' '
      ELSE
         CALL file (opfile, klist, 'lls')
         WRITE(*,1050)klist

         WRITE(*,*)'     CALLING read_klbllist'
         CALL read_klbllist (klist, ngas, idgas, isogas, nwave, vwave)
         WRITE(*,*)'     read_klist COMPLETE'
         WRITE(*,*)' '
      ENDIF


C Now read the scattering files if required.
      scatterf = .FALSE.
      DO I=1,npath
        IF(imod(I).EQ.15.OR.imod(I).EQ.16)scatterf = .TRUE.
      ENDDO

      IF(scatterf)THEN
        CALL file(runname, radfile, 'sca')
        CALL get_scatter(radfile,ncont)
      ENDIF

C ... and the xsc files likewise.
      dustf = .FALSE.
      IF(ncont.GT.0)dustf = .TRUE.

      IF(dustf)THEN
        CALL file(runname, xscfil, 'xsc')
        CALL get_xsec(xscfil, ncont)
      ENDIF

      PRINT*,'NPATH = ',NPATH

C=======================================================================
C
C	Call CIRSradg_wave.	
C
C=======================================================================

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
      WRITE(*,*)' '
      WRITE(*,*)' Waiting for CIRSradg_wave etal to complete ...'
      xcomp = 0.0
      DO 1000 iwave=1,nwave
        if(nwave.gt.1)then
          xref = 100.0*FLOAT(iwave - 1)/FLOAT(nwave - 1)
        else
          xref = 0.0
        endif
        IF(xref.GE.xcomp)THEN
          WRITE(*,*)' Percent Complete = ',xcomp
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
    

	CALL cirsradg_wave (dist,inormal,iray,delh,nlayer,npath,ngas,
     1  press,temp,pp,amount,iwave,ispace,AMFORM,vv,nlayin,
     2  layinc,cont,scale,imod,idgas,isogas,emtemp,itype1,
     3  nem, vem, emissivity, tsurf, gtsurf, RADIUS1,
     4  flagh2p,hfp,radextra,nsw,isw,output,doutputdq)


C Convert from rates of change with respect to layer variables to rates
C of change of .prf profile variables.
        CALL map2pro(nparam,doutputdq,doutmoddq)


C Convert from rates of change of .prf profile variables to rates of
C change of desired variables.
        CALL map2xvec(nparam,npro,npath,nv,xmap,doutmoddq,doutdx)

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
1000  CONTINUE

C Convolve output calculation spectra with a square bin of width FWHM to
C get convoluted spectra.
      DO ipath=1,npath

        DO iwave=1,nwave
           ioff1 = nwave*(ipath - 1) + iwave
           y(iwave) = tempout(ioff1)
           ygt(iwave) = tempgtsurf(ioff1)
        ENDDO		

        CALL cirsconv(runname,fwhm1,nwave,vwave,y,nconv,vconv,yout)	
        CALL cirsconv(runname,fwhm1,nwave,vwave,ygt,nconv,vconv,
     1			youtgt)	
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

          CALL cirsconv(runname,fwhm1,nwave,vwave,y,nconv,vconv,yout)

          DO iconv=1,nconv
            ioff2 = nconv*nv*(ipath - 1) + (iv - 1)*nconv + iconv
            gradients(ioff2) = yout(iconv)
C            print*,'YYY',iconv,iv,ipath,ioff2,gradients(ioff2)
          ENDDO
        ENDDO
      ENDDO

C=======================================================================
C
C	Wrap up: formats, return, and end.
C
C=======================================================================

      CALL close_scat(ncont)

1050  FORMAT (/'Klist filename: ', A)

      RETURN

      END
