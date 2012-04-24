      subroutine get_hg(vwave,calpha,ncont,icont,dphase)
C***********************************************************************
C_TITL:	GET_HG
C
C_DESC:	
C
C_ARGS:	Input variables:
C	vwave		REAL		Wavenumber.
C	calpha		REAL		Cos(scattering angle)
C	ncont		INTEGER		Number of dust types included.
C	icont		INTEGER		Number of dust types index.
C
C	Output variables:
C	dphase		DOUBLE		Resulting phase function
C
C_CALL: No calls.
C
C_HIST:	PGJI	Original version.
C***********************************************************************

      IMPLICIT NONE
      include '../includes/arrdef.f'
      REAL vwave
      REAL xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec),
     1 xg2(maxcon,maxsec)
      REAL twave,tico,x,y,z,frac,xw,tnco
      DOUBLE PRECISION dphase,calpha,f,g1,g2,a,b
      INTEGER iunit,ico,i,ncont,icont,nco,j
C IUNIT: File unit
      CHARACTER*100 ipfile

      COMMON /hgphas/ xwave,xf,xg1,xg2,tnco,twave,frac,tico

C********************************* CODE ********************************

      iunit = 11

      IF(xwave(1).LT.0.0)THEN
        DO 50 i=1,ncont
          ipfile = 'hgphase*.dat'
          ipfile(8:8) = char(i+48)
cc          WRITE(*,*)' GET_HG.f :: Reading ipfile: ',ipfile
          ico = 0
          OPEN(iunit,FILE=ipfile,STATUS='old')
10        READ(iunit,*,END=99)xw,x,y,z
          ico=ico+1
          if(ico.gt.maxsec)then
           print*,'Error in get_hg: too many wavelengths',maxsec
           stop
          endif
          xwave(ico) = xw
          xf(i,ico) = x
          xg1(i,ico) = y
          xg2(i,ico) = z
          GOTO 10
99        CONTINUE
          CLOSE(iunit)
          nco = ico
cc          WRITE(*,*)' GET_HG.f :: nco = ',nco
          tnco = nco+0.1
50      CONTINUE
      ELSE
        nco = INT(tnco)
      ENDIF

      IF(vwave.NE.twave)THEN
        ico = 0
        DO 60 i=1,nco-1
          IF(xwave(i).LE.vwave.AND.xwave(i+1).GT.vwave)ico = i
60      CONTINUE
        IF(vwave.EQ.xwave(nco))ico = nco-1
        IF(ico.EQ.0)THEN
          WRITE(*,*)' GET_HG.f :: Error: wavenumber out of range.'
          WRITE(*,*)' GET_HG.f :: Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' GET_HG.f :: vwave: ',vwave
          WRITE(*,*)' GET_HG.f :: icont,ncont: ',icont,ncont
          WRITE(*,*)' GET_HG.f :: xwave ,xf, xg1, xg2: '
          DO i=1,nco
            WRITE(*,*)xwave(i),xf(icont,i),xg1(icont,i),xg2(icont,i)
          ENDDO
          STOP
        ENDIF
        frac = (vwave - xwave(ico))/(xwave(ico+1)-xwave(ico))
        tico = ico + 0.1
cc        WRITE(*,*)'vwave,xwave(ico),xwave(ico+1),frac,ico'
cc        WRITE(*,*)vwave,xwave(ico),xwave(ico+1),frac,ico
        twave = vwave
      ELSE
        ico = int(tico)
      ENDIF
  
      f = (1.0 - frac)*xf(icont,ico)  + frac*xf(icont,ico+1)
      g1 = (1.0 - frac)*xg1(icont,ico) + frac*xg1(icont,ico+1)
      g2 = (1.0 - frac)*xg2(icont,ico) + frac*xg2(icont,ico+1)

      a = (1 - g1**2)/(1 + g1**2 - 2*g1*calpha)**1.5 
      b = (1 - g2**2)/(1 + g2**2 - 2*g2*calpha)**1.5 

      dphase = f*a + (1 - f)*b

      RETURN

      END
