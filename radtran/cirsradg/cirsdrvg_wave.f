	PROGRAM CIRSdrvg_wave
C*********************************************************************
C_TITL:	CIRSdrvg_wave
C
C_DESC:	Driver program for new gradient-version of cirsrad.
C
C_ARGS:	
C
C_CALL: gettime		Obtains the system time for use in determining the
C                       total time used by the program since execution.
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

      INTEGER intmod,nconv,nwave,i,j,itype
      INTEGER iconv,ioff1,ioff2,iv,nv,iwave,ipath,ispace
      INTEGER npoints,planet,inormal,iray,iptf,npath,ichannel,nem

      REAL vconv(maxbin),vwave(maxbin)
      REAL vref,delvk,vmin,vmax,vem(MAXSEC),emissivity(MAXSEC),tsurf
      REAL dist,wnumbot,wnumtop,delv,fwhm,gtsurf

      REAL tot_time
      DOUBLE PRECISION time,time1,time2

      REAL calcout(maxout3),gradients(maxout4)
      REAL xmap(maxv,maxgas+2+maxcon,maxpro)
      REAL y(maxout)
      CHARACTER*80 text
      CHARACTER*100 runname,patfile,outfile,ciafil
      CHARACTER*1 ANS

      REAL xwave(100),xf(10,100),xg1(10,100),xg2(10,100)
      REAL tnco,twave,frac,tico
      COMMON /hgphas/ xwave,xf,xg1,xg2,tnco,twave,frac,tico

      INCLUDE '../includes/ciacom.f'
      INCLUDE '../includes/gascom.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA


C-----------------------------------------------------------------------
C
C	Write a header
C
C-----------------------------------------------------------------------

C Obtain the system time for use in determining the total time used by the
C program since execution.
      CALL gettime(time)
      time1= time

C Reset xwave to force read of hgphase* files for scloud10-11.
      xwave(1)=-1

      WRITE(*,*)' '
      WRITE(*,*)'           WELCOME TO CIRSDRVG_WAVE'
      WRITE(*,*)' '

C-----------------------------------------------------------------------
C
C	Prompt the user for prescribed variables.
C
C-----------------------------------------------------------------------

      WRITE(*,*)'Give operation filename: '
      READ(*,1000)runname
1     FORMAT(A)
      CALL FILE(runname,patfile,'pat')
      OPEN(2,FILE=PATFILE,STATUS='OLD')
2     READ(2,1,END=9)TEXT
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
      IF(TEXT(1:1).EQ.' ')THEN	! skip blank line
        GOTO 2
      ELSE IF(TEXT(1:8).EQ.'INTERVAL')THEN
        READ(2,*)wnumbot,wnumtop,delv,fwhm
      ENDIF
      GOTO 2

9     CONTINUE

      CLOSE(2)

C Calculate number of points in output array 
      npoints= NINT((wnumtop-wnumbot)/delv) + 1

      CALL PROMPT('Enter ISPACE : ')
      READ*,ispace

      nconv= npoints
      WRITE(*,*)' CIRSDRVG_WAVE.f :: # of convolution wavelengths: ',
     1 nconv
      WRITE(*,*)' CIRSDRVG_WAVE.f :: i, convolution wavelengths: '
      DO i= 1, nconv
        vconv(i)= wnumbot + ((i-1) * DELV)
        WRITE(*,*)i,vconv(i)
      ENDDO

      print*,'OK (Y/N)?'
      read(5,1)ans
      ichannel=0
      if(ans.ne.'y'.and.ans.ne.'Y')then
       call prompt('Enter nconv : ')
       read*,nconv
       print*,'Enter wavenumbers:'
       read*,(vconv(i),i=1,nconv)
       ichannel=1
       fwhm=0.0
      endif


10    WRITE(*,*)'Select a planet (1-9) to look up solar dist: '
      READ (*,*) Planet
      IF (planet.EQ.1) Dist= 0.39	! Mercury
      IF (planet.EQ.2) Dist= 0.72	! Venus
      IF (planet.EQ.3) Dist= 1.00	! Earth
      IF (planet.EQ.4) Dist= 1.50	! Mars
      IF (planet.EQ.5) Dist= 5.20	! Jupiter
      IF (planet.EQ.6) Dist= 9.50	! Saturn
      IF (planet.EQ.7) Dist= 19.2	! Uranus
      IF (planet.EQ.8) Dist= 30.1	! Neptune
      IF (planet.EQ.9) Dist= 39.5	! Pluto
      IF ((planet.GT.9).OR.(planet.EQ.0)) GOTO 10


C-----------------------------------------------------------------------
C
C       Two wavenumber arrays are considered; one for final output
C       (VConv), and one for actual use in calculating outputs (VWave).
C
C	NCONV	the number of data points for the prompted wavenumber
C		range and resolution.
C	VCONV	an array containing each wavenumber value for use when
C		convolving.
C
C	NWAVE	the number of points in the k-table that cover the
C		prompted wavenumber range.
C	VWAVE	an array of the points within this wavenumber range.
C
C	There are two ways of producing NWAVE & VWAVE: fast and slow.
C	In the slow mode, NWAVE & VWAVE are the same as the convolution
C	wavenumbers. In the fast mode, they are the same wavenumbers as
C	the K tables are calculated at.
C
C	Accuracy of methods is approximately the same.
C
C-----------------------------------------------------------------------

      if(ichannel.eq.1)then
       intmod=1
      else
30     WRITE(*,*)'Enter Wavenumber interpolation mode:'
       WRITE(*,*)'    (1) Fast post-calculation'
       WRITE(*,*)'    (2) Slow K-coeff regridding'
       WRITE(*,*)'Option: '
       READ(*,*)intmod
       IF ((intmod.NE.1).AND.(intmod.NE.2)) GOTO 30
      endif
      WRITE(*,*)' '
      IF (intmod.EQ.1) WRITE(*,*)'-> Fast post-calculation selected'
      IF (intmod.EQ.2) WRITE(*,*)'-> Slow K-coeff regridding selected'
      WRITE(*,*)' '

      IF (intmod.EQ.1) THEN
        if(ichannel.eq.0)then
         WRITE(*,*)'Enter reference wavenumber and step in ktables'
         READ (*,*) VREF,DelVK

         vmin = vref + delvk*INT((wnumbot - 0.5*fwhm-vref)/delvk)
         vmax = vref + delvk*(NINT((wnumtop + 0.5*fwhm-vref)/delvk) + 1)

C Calculate number of points in tabulated wavenumber array plus one to
C account for the end point
         nwave= NINT((vmax - vmin) / delvk) + 1
         DO j= 1, nwave
           vwave(j)= VMIN + ((j-1) * delvk)
         ENDDO
        else
         nwave=nconv
         do j=1,nwave
          vwave(j)=vconv(j)
         enddo
        endif
      ELSE
	nwave= nconv + 2
	vwave(1)= vconv(1) - delv
	vwave(nwave)= vconv(nconv) + delv
	DO i= 1, nconv
	  vwave(i+1)= vconv(i)
	ENDDO
      ENDIF

      WRITE(*,*)' CIRSDRVG_WAVE.f :: # of calculation wavelengths: ',
     1  NWave
      WRITE(*,*)' CIRSDRVG_WAVE.f :: i, calculation wavelengths: '
      DO i=1,nwave
        WRITE(*,*)i,vwave(I)
      ENDDO

C-----------------------------------------------------------------------
C
C       Get scattering type and Ortho/Para information.
C
C-----------------------------------------------------------------------
      itype=11


      CALL READFLAGS(RUNNAME,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF)

      CALL FILE(runname,CIAFIL,'cia')

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
	
      nv=1
      DO i=1,300
        xmap(1,12,i)=1.0
      ENDDO

      CALL cirsrtfg_wave (runname, dist, inormal, iray, fwhm, ispace, 
     1 vwave,nwave,itype, nem, vem, emissivity, tsurf, gtsurf, nv, xmap, 
     2 vconv, nconv, npath, calcout, gradients)

      CALL FILE(runname,outfile,'out')
      OPEN(2,FILE=outfile,STATUS='unknown')
      WRITE(2,*)npath,nconv,nv
      DO ipath=1,npath
	WRITE(2,*)ipath
        DO iconv=1,nconv
          ioff1 = nconv*(ipath-1) + iconv
          WRITE(2,*)vconv(iconv),calcout(ioff1)
        ENDDO
        DO iv=1,nv
          DO iconv=1,nconv
            ioff2 = nconv*nv*(ipath - 1) + (iv - 1)*nconv + iconv
            y(iconv)=gradients(ioff2)
          ENDDO
          WRITE(2,*)(y(iconv),iconv=1,nconv)
        ENDDO
      ENDDO

C----------------------------------------------------------------------- 
C 
C	Wrap up: formats and end.
C
C-----------------------------------------------------------------------

      WRITE(*,*)' CIRSDRVG_WAVE.f :: calculation complete.'
      CALL GETTIME(time)
      time2 = time
      tot_time = sngl(time2-time1)
      WRITE(*,200)tot_time
200   FORMAT(' Elapsed time including convolution (sec)= ',F8.1)

1000  FORMAT (A60)

      END
