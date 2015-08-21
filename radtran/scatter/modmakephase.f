      subroutine modmakephase(iwave,imode,inorm,iscat,
     1   parm,rs,srefind,runname,lambda0)
C     ****************************************************************
C     Subroutine to allow Makephase functionality from within other programs
C
C     Input variables
C	iwave	integer	Defines wavelength space. 1=wavelength, 2=wavenumber
C	imode	integer	ID of particle type to adjust
C	inorm	integer	Flag to see if phase function is to be renormalise.
C			   1 = renormalise
C	iscat	 integer	Scattering type
C	parm(3) 	real	Scattering parameters
C	rs(3)		real	Size distribution of particles
C	srefind(max_wave,2) real complex refractive index spectra
C	outfile	character*100	Name of associated .xsc file
C	lambda0	real	normalising extinction x-section wavelength 
C				(-ve to disable)
C     Output variables
C	None
C
C     Modified from Makephase. Pat Irwin.  20/2/14
C
C     ****************************************************************
C
	implicit none
	integer	max_theta, max_mode, max_wave,inorm,imode
	parameter (max_theta = 100)
	parameter (max_mode = 10)
	parameter (max_wave = 1000)

	integer	nmode, iwave, nwave, ntheta, nphase, I, iscat, 
     1		irec, J, siscat(max_mode), wc, K, jref
	real	pi, wmin, wmax, wstep, minlam, lambda, parm(3), 
     1		rs(3), refind(2), w,  theta(max_theta), scat, 
     2		ext, phase(max_theta), omega, calpha,  x, qext, gsec,
     3		xsec(max_mode,max_wave,2), qsca, qabs, sum, du, dsum,
     5		scalef,f,g1,g2,rms,srefind(max_wave,2),
     6          phase1(max_theta),wave(max_wave),lambda0,xx,yy
	double precision f1, f2, hg11, hg12, hg21, hg22
	complex	nc

	character*1 	ans
	character*10	wavetype(2)
	character*100 	runname,outfile,hgfile
	character*512 	buffer
      

        integer         ilist, nspec, idspec
        real          table(5000,3)
        common /store/ ilist,nspec,idspec,table

C       Extra variables to make sure that hgphase files are read in again
C       by read_hg.f
        integer mcon,mphas
        parameter(mcon=10,mphas=300)
        real xwave(mphas),xf(mcon,mphas),xg1(mcon,mphas)
        real xg2(mcon,mphas),tnco,twave,frac,tico
        common /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico


	pi = 3.141592654
	wavetype(1) = 'wavelength'
	wavetype(2) = 'wavenumber'

        print*,'modmakephase'
        print*,'iwave,imode,inorm,iscat,lambda0',iwave,imode,
     1  inorm,iscat,lambda0
        print*,'(parm(i),i=1,3)',(parm(i),i=1,3)
        print*,'(rs(i),i=1,3)',(rs(i),i=1,3)

        call get_xsecA(runname,nmode,nwave,wave,xsec)

        print*,'Refindex'
        print*,'nwave = ',nwave
        do i=1,nwave
         print*,wave(i),srefind(i,1),srefind(i,2)
        enddo

	call file (runname,runname,'xsc')
        open (12, file=runname, status='unknown')

	call get_theta (theta, max_theta, ntheta, nphase)

	if (iwave.eq.1) then
		minlam = wave(1)
	else
		minlam = 1.e4/wave(nwave)	! From wavenumbers
	endif

        if(imode.lt.0.or.imode.gt.nmode)then
         print*,'Error in modmakephase.f. imode is out of range'
         print*,'imode,nmode',imode,nmode
        endif

C-----------------------------------------------------------------------
C
C	Begin calculation over number of modes. Get relevant inputs. 
C	Write mode dependent headers to the PHASEN files (N is replaced 
C	by the number of the mode under consideration).
C
C-----------------------------------------------------------------------

 	ilist=0

C                print*,iscat,parm,rs

	outfile = 'PHASEN.DAT'
	outfile(6:6) = char(imode+48)
	open (13,file=outfile, status='unknown', 
     1        access='direct',recl=512, form='formatted')
	irec = 1
        wmin=wave(1)
        wmax=wave(nwave)
        wstep=wave(2)-wave(1)
	write (buffer, 1010) wavetype(iwave), wmin, wmax, wstep,
     1				nwave, nphase
	write (13, 1000, rec=irec) buffer
	irec = irec + 1

	if (iscat.eq.1) then
		write (buffer, 1020) (parm(J), J = 1, 3),
     1			(rs(J), J = 1, 3)
	elseif (iscat.eq.2) then
		write (buffer, 1030) (parm(J), J = 1, 2),
     1			(rs(J), J = 1, 3)
	elseif (iscat.eq.3) then
		write (buffer, 1035) (parm(J), J = 1, 3),
     1			(rs(J), J = 1, 3)
	elseif (iscat.eq.4) then
		write (buffer, 1040) parm(1)
	elseif (iscat.eq.5) then
		write (buffer, 1050) 
	elseif (iscat.eq.6) then
		write (buffer, 1060) (parm(J), J = 1, 3)
	elseif (iscat.eq.7) then
		write (buffer, 1070) parm(1)
	endif

	write(13,1000,rec=irec) buffer
	irec = irec + 1

	write (buffer, 1080) (theta(J), J = 1, nphase)
	write(13,1000,rec=irec) buffer
	irec = irec + 1


	hgfile = 'hgphaseN.dat'
	hgfile(8:8) = char(imode+48)
	open (14,file=hgfile, status='unknown')
                

C-----------------------------------------------------------------------
C
C	Begin loop over wavenumbers or wavelengths.
C
C-----------------------------------------------------------------------

	do J = 1, nwave
		w = wave(J)
C                       print*,i,w
		if (iwave.eq.1) then
			lambda = w
		else
			lambda = 1.e4/w
		endif

                refind(1)=srefind(J,1)
                refind(2)=srefind(J,2)

C-----------------------------------------------------------------------
C
C	Do Mie scattering first. Omega is the single scattering albedo.
C
C-----------------------------------------------------------------------

		if (iscat.lt.5) then
			call miescat(lambda, iscat,parm, 1, rs,
     1				refind, theta, 
     2				ntheta, scat, ext, phase, 
     3			        nphase)
			omega = scat/ext

		end if


C-----------------------------------------------------------------------
C
C	Do isotropic scattering.
C
C-----------------------------------------------------------------------

		if (iscat.eq.5) then	
			do K = 1, nphase
				phase(K) = 1.0
         		end do
		endif
		if ((iscat.eq.5).or.(iscat.eq.6)) then
			if (iwave.eq.1) then
				write (*,* ) 
     1					' Wavelength is: ', w
			else
				write (*,* ) 
     1					' Wavenumber is: ', w
			endif
  			write (*,'('' Give extinction'',
     1  '' cross-section (cm2) and single scattering albedo: '',$)')
				read (*,*) ext, omega
		endif

C-----------------------------------------------------------------------
C
C	Do Henyey-Greenstein scattering.
C
C-----------------------------------------------------------------------

		if (iscat.eq.6) then
			f1 = parm(1)
			f2 = 1.0d0-f1
			hg11 = 1.0d0-parm(2)*parm(2)
			hg12 = 2.0d0-hg11
			hg21 = 1.0d0-parm(3)*parm(3)
			hg22 = 2.0d0-hg21
			do K = 1, nphase	
				calpha=cos(theta(K)*pi/180)
				phase(K) = sngl(f1 * hg11 / 
     1	dsqrt( hg12 - 2.0d0*parm(2)*calpha )**3 + f2 * hg21 / 
     2	dsqrt( hg22 - 2.0d0*parm(3)*calpha )**3)
				end do
		endif

C-----------------------------------------------------------------------
C
C	Do dipole scattering.
C
C-----------------------------------------------------------------------

		if (iscat.eq.7) then	
			nc=cmplx(refind(1),-refind(2))
			x = 2 * pi * parm(1)/lambda
			qsca = (8./3.)*(x**4) * abs((nc**2-1)/
     1				(nc**2 + 2))
			qabs = - 4 * x * aimag((nc**2 - 1)/
     1				(nc**2 + 2))
			qext = qsca + qabs
			omega = qsca/qext
			gsec = pi * (parm(1) * 1.e-4)**2
			ext = qext * gsec
   
			do K = 1, nphase
				calpha = cos(theta(K)*pi/180)
 				phase(K) = 0.75 * (1.0 + 
     1					calpha*calpha)
			end do
		end if 
  
C-----------------------------------------------------------------------
C
C	Write output data at each wavenumber or wavelength, then close 
C	PHASEN file. Also perform renormalization to ensure that the
C       intgeral of the phase function from cos(theta) = -1 to +1 is equal
C       to 1/2pi
C
C-----------------------------------------------------------------------

	 do K = 1, nphase
		phase1(K)=phase(k)
	 	phase(K) = 0.25 * phase(K)/pi
	 enddo


         sum=0
         do k=1,nphase-1
              du = cos(theta(k)*pi/180.0) - 
     &                        cos(theta(k+1)*pi/180.0)
                   sum=sum+0.5*(phase(k)+phase(k+1))*du
              end do
              sum=sum*2*pi
C             print*,'modmakephase: 2*pi*Int(pdu) = ',sum
              dsum = abs(sum - 1.0)
              if(dsum.gt.0.05)then
                     print*,'modmakephase: WARNING, > 5% off unity'       
              end if
              if(inorm.eq.1)then
C                   print*,'Renormalising : '
                     do k=1,nphase
                        phase(k)=phase(k)*1.0/sum
                        phase1(k)=phase1(k)*1.0/sum
                     end do
              endif

              if(phase(1).ge.10000.0)then
 	        write (buffer,1095,ERR=666) 
     &		  w, ext, omega,(phase(K),K=1, nphase)
	      else
		write (buffer,1090,ERR=666) 
     &		   w, ext, omega,(phase(K),K=1, nphase)
	      endif
              write(13,1000,rec=irec)buffer

	      goto 667
 666	      print *, 'WARNING: error occured during write',
     &                   ' to buffer. Probably the forward'
	      print *, 'scattering amplitude exceeded the',
     &                  ' FORMAT limits due to large particle size.'

 667	      irec=irec+1

C	      Subfithgm is expecting phase functions normalised
C	      to 4pi, rather than 1, so pass phase1 instead of
C	      phase.
	      call subfithgm(nphase,theta,phase1,f,g1,g2,rms)



              write(14,*)w,f,g1,g2

	      xsec(imode,J,1) = ext
	      xsec(imode,J,2) = omega

	enddo

C       Since the hgphase file has bee updated we need to make sure that
C       read_hg.f reads it in again, by setting xwave(1) < 0.
        xwave(1)=-1.


 	close(13)
 	close(14)

C-----------------------------------------------------------------------
C
C	Now produce .XSC files
C
C-----------------------------------------------------------------------

        if(lambda0.gt.0)then
         xx=1e10
         jref=1
         do i=1,nwave
          yy=abs(wave(i)-lambda0)
          if(yy.lt.xx)then
           xx=yy
           jref=i
          endif
         enddo
         do j=1,nmode
          xx=xsec(j,jref,1)
          do i=1,nwave
           xsec(j,i,1)=xsec(j,i,1)/xx
          enddo
         enddo
        endif

	write (12,*) nmode        
	do I = 1, nwave
		write(12,*) wave(i),(xsec(J,I,1), J=1, nmode)
		write(12,*)(xsec(J,I,2), J=1, nmode)
	enddo

        close(12)

C-----------------------------------------------------------------------
C
C	Format statements and end.
C
C-----------------------------------------------------------------------

1000	format (a512)
1010	format (1x, a10, 2(2x, f8.2), 2x, f8.4, 2(2x, i4))
1020	format (' Mie, standard gamma, A = ',0pf7.2, ' B = ',0pf7.2,
     1   ' alpha = ',0pf7.2, ' rmin = ', 0pf7.3, ' rmax = ',
     2   0pf7.2, ' rstep = ', 0pf7.3)
1035	format (' Mie, MCS gamma, A = ',0pf7.2, ' B = ',0pf7.2,
     1   ' C = ',0pf7.2, ' rmin = ', 0pf7.3, ' rmax = ',
     2   0pf7.2, ' rstep = ', 0pf7.3)
1030	format (' Mie, log normal ', ' R0 = ',0pf7.2, ' Sig = ',0pf7.2,
     1   ' rmin = ', 0pf7.3, ' rmax = ',0pf7.2, ' rstep = ', 
     2   0pf7.3)	
1040	format (' Mie, single particle ', ' R0 = ',0pf7.2)	
1050	format (' Isotropic scattering')
1060	format (' Henyey-Greenstein, f = ', 0pf7.2, ' g1 = ',0pf7.2,
     1   ' g2 = ',0pf7.2)
1070	format (' Dipole, R0 = ', 0pf7.2)
1080	format (100(1x,f8.3))
1090    format (1x,f8.2,1x,e12.5,1X,e12.5,50(1x,f10.5))
1095    format (1x,f8.2,1x,e12.5,1X,e12.5,50(1x,f10.4))

	end


