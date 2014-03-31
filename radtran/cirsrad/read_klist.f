C***********************************************************************
C_TITL:	READ_KLIST
C
C_DESC:	Writes variables to common block INTERPK for subsequent use in 
C	CIRSRADG. Variables passed from the main routine are:
C
C_ARGS:	Input variables:
C	klist		CHARACTER*100	.kls file containing names of
C					.kta files to be used.
C	ngas		INTEGER		Number of gases
C	idgas(ngas)	INTEGER		Gas IDs
C	isogas(ngas)	INTEGER		Isotope IDs
C	nwave		INTEGER		Number of wavenumbers
C	vwave		REAL		Wavenumbers for calculation
C
C_FILE:	unit=11		klist (<runname>.kls)
C
C_CALL:	FILE		Forces file extension.
C	READ_KHEAD	Opens .kta file, reads in header and returns
C			the beginning record number for the user-selected
C			spectral range as well as temp and pressure ranges
C			used.
C	FINDLOC		For a given monotonic array and 2 values, finds
C			the indices of the array values contained by the
C			passed variables.
C
C_HIST:	???????	???	Orginal version.
C***********************************************************************

      SUBROUTINE read_klist (klist,ngas,idgas,isogas,nwave,vwave)

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
C Defines the maximum values for a series of variables (layers, bins,
C paths, etc.)

      INTEGER maxkfil,imatch,ipo1
      real MAXDX,MAXDX1
      PARAMETER (maxkfil=100,MAXDX=0.01)

      INTEGER ngas,nwave,nkl,npk,ntk,ngk,i,j,i1,i2,k,npoint
C NPK: Number of pressure points in k-tables.
C NTK: Number of temperature points in k-tables.
C NGK: Number of g-ordinates in k-tables.
      INTEGER idgas(ngas),isogas(ngas)
      INTEGER unit(maxkfil),npt(maxkfil),idgask(maxkfil),
     1 isogask(maxkfil),np(maxkfil),nt(maxkfil),ng(maxkfil),
     2 irec(maxkfil)
      INTEGER iflag(maxbin,maxgas),lun(maxbin,maxgas),
C LUN: Unit number of file. Set to -99 for missing points.
     1 ireck(maxbin,maxgas),ipo
C IRECK: Beginning record number in file.

      REAL pmin,pmax,tmin,tmax,fwhmk,delvk
      REAL vwave(nwave),vcen(mpoint),xcenk(mpoint)
      REAL xmin(maxkfil),xmax(maxkfil),delx(maxkfil),fwhm1(maxkfil)
      REAL pk(maxk),tk(maxk)
C PK: Pressure values in k-tables.
C TK@ Temperature values in k-tables.
      REAL g_ord(maxg),delg(maxg)
C G_ORD: g-ordinates.
C DELG: g-weights.
      REAL xmink(maxbin,maxgas),delk(maxbin,maxgas),DX
C XMINK: Beginning wavenumber in table.
C DELK: Wavenumber step of evaluation points.
      REAL kout(maxlay,maxgas,maxg),dkoutdt(maxlay,maxgas,maxg)
      CHARACTER*100 klist,ktafil(maxkfil),null

      COMMON /interpk/ lun,ireck,xmink,delk,pk,npk,tk,ntk,ngk,
     1 delvk,fwhmk,g_ord,delg,kout,dkoutdt

C********************************* CODE ********************************

C-----------------------------------------------------------------------
C
C	Read in names of files.
C
C-----------------------------------------------------------------------

      nkl = 0
      null = ' '
      print*,'klist file : ',klist
      OPEN (UNIT=11,FILE=klist,STATUS='old')

10    nkl = nkl + 1
      READ(11,1020,END=20)ktafil(nkl)
      IF(ktafil(nkl).EQ.null)GOTO 20

      CALL file(ktafil(nkl),ktafil(nkl),'kta')
      WRITE(*,1030)ktafil(nkl)
      GOTO 10
	
20    nkl = nkl - 1
      CLOSE(UNIT=11)

C-----------------------------------------------------------------------
C
C	Open appropriate Correlated K files and compare. 
C
C-----------------------------------------------------------------------
      WRITE(*,*)' READ_KLIST.f :: Number of k-files = ',nkl

      DO i=1,nkl
        unit(i) = 100 + i
        WRITE(*,1035)KTAFIL(i)
        CALL read_khead(ktafil(i),unit(i),npt(i),xmin(i),delx(i),
     1  fwhm1(i),vcen,idgask(i),isogask(i),pk,tk,np(i),nt(i),g_ord,
     2  delg,ng(i),irec(i))
        xmax(i) = xmin(i) + (npt(i) - 1) * delx(i)

        IF (i.EQ.1) THEN
          pmin = pk(1)
          pmax = pk(np(i))
          tmin = tk(1)
          tmax = tk(nt(i))
          npoint = npt(i)
          if(npoint.gt.mpoint)then
           print*,'Error in READ_KLIST: npoint > mpoint'
           print*,npoint,mpoint
           print*,'Use smaller k-tables or modify MPOINT in arrdef.f'
           print*,'and recompile everything'
           stop
          endif
          npk = np(i)
          ntk = nt(i)
          ngk = ng(i)
C          print*,'read_klist : ngk = ',ngk
          delvk = delx(i)
          fwhmk = fwhm1(i)
          do j=1,npoint
           xcenk(j)=vcen(j)
          enddo
        ELSE
          IF((np(i).NE.npk).OR.(nt(i).NE.ntk).OR.(ng(i).NE.ngk)
     1    .OR.(fwhm1(i).NE.fwhmk).OR.(delx(i).NE.delvk)) THEN
            WRITE(*,*)' READ_KLIST.f :: Error: Problems reading'
            WRITE(*,*)' k-tables array sizes. Stopping program.'
            WRITE(*,*)' '
            WRITE(*,*)' ktafil(i) ',ktafil(i)
            WRITE(*,*)' npk = ',npk,' np(i) = ',np(i)
            WRITE(*,*)' ntk = ',ntk,' nt(i) = ',nt(i)
            WRITE(*,*)' ngk = ',ngk,' ng(i) = ',ng(i)
            WRITE(*,*)' fwhmk = ',fwhmk,' fwhm(i) = ',fwhm1(i)
            WRITE(*,*)' delvk = ',delvk,' delx(i) = ',delx(i)
            STOP
          ENDIF

          IF(delx(i).lt.0)THEN
           IF(npt(i).NE.npoint)THEN
            WRITE(*,*)'Error in READ_KLIST: npoints dont match'
            WRITE(*,*)' npoint = ',npoint,' npt(i) = ',npt(i)
            STOP
           ENDIF
           DO j=1,npoint
            dx = 100*ABS((VCEN(J)-XCENK(J))/VCEN(J))
            IF(dx.gt.MAXDX)THEN
             WRITE(*,*)'Error in READ_KLIST: wave centres dont match'
             WRITE(*,*)J,VCEN(J),XCENK(J)
             STOP
            ENDIF
           ENDDO
          ENDIF

          IF ((abs(pk(1)-pmin).GT.1.e-4).OR.
     1    (abs(pk(npk)-pmax).GT.1.e-4).OR.
     2    (abs(tk(1)-tmin).GT.1.e-4).or.
     3    (abs(tk(ntk)-tmax).GT.1.e-4)) THEN
            WRITE(*,*)' READ_KLIST.f :: Error: Problems reading'
            WRITE(*,*)' k-tables PT Grid. Stopping program.'
            WRITE(*,*)' '
            WRITE(*,*)' ktafil(I) ',ktafil(I)
            WRITE(*,*)' pmin = ',pmin,' pk(1) = ',pk(1)
            WRITE(*,*)' pmax = ',pmax,' pk(npk) = ',pk(npk)
            WRITE(*,*)' tmin = ',tmin,' tk(1) = ',tk(1)
            WRITE(*,*)' tmax = ',tmax,' tk(ntk) = ',tk(ntk)
            STOP
          ENDIF

        ENDIF
      ENDDO
	
C-----------------------------------------------------------------------
C
C	Check some variable sizes.
C
C-----------------------------------------------------------------------
C      print*,'test'

      IF (ngk.GT.maxg) THEN
        WRITE(*,*)' READ_KLIST.f :: Error: Too many g-ordinates.'
        WRITE(*,*)' Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)' ngk = ',ngk,' maxg = ',maxg
        STOP
      ENDIF

      IF ((npk.GT.maxk).OR.(ntk.GT.maxk)) THEN
        WRITE(*,*)' READ_KLIST.f :: Error: Too many P/T points in'
        WRITE(*,*)' k-tables. Stopping program.' 
        WRITE(*,*)' '
        WRITE(*,*)' maxk = ',maxk
        WRITE(*,*)' np = ',np,' nt = ',nt
        STOP
      ENDIF

C-----------------------------------------------------------------------
C
C	Assign relevant correlated K files by wavenumber and gas. 
C	For H2, He, H2S assign dummy values. H2 and He are calculated as
C	continuum contributions only. Iflag is used so that variables set
C	by the first K table encountered are not reset, (thus overlapping
C	K tables can be used, their priorities decided by the order in
C	which they appear in the KLS file.)
C
C-----------------------------------------------------------------------

      DO i=1,ngas
        DO j=1,nwave
          iflag(j,i) = 0
        ENDDO
      ENDDO

      DO i=1,ngas
        IF ((idgas(i).EQ.39).OR.(idgas(i).EQ.40).OR.
     1	(idgas(i).EQ.22).OR.(idgas(i).EQ.36)) THEN
          DO j = 1, nwave
            iflag(j,i) = 1
            lun(j,i) = 0
            xmink(j,i) = 0.0
            delk(j,i) = 0.0
            ireck(j,i) = 0
          ENDDO
          GOTO 30
        ENDIF
        DO j=1,nkl
          IF ((idgas(i).EQ.idgask(j)).AND.
     1	  (isogas(i).EQ.isogask(j))) THEN
            IF(delx(j).gt.0.0)THEN
              CALL findloc(vwave,nwave,xmin(j),xmax(j),i1,i2)
              IF((i1.NE.0).AND.(i2.NE.0))THEN
                DO k=i1,i2
                  IF (iflag(k,i).EQ.0) THEN
                    iflag(k,i) = 1
                    lun(k,i) = unit(j)
                    xmink(k,i) = xmin(j)
                    delk(k,i) = delx(j)
                    ireck(k,i) = irec(j)
                  ENDIF
                ENDDO
              ENDIF
            ELSE
             DO k=1,nwave
              imatch=0
              DO ipo=1,npoint
               dx = 100*ABS((vwave(k)-XCENK(ipo))/XCENK(ipo))
               IF(dx.lt.MAXDX)THEN
C                  print*,'read_klist',k,vwave(k),ipo,xcenk(ipo)
                  imatch=1
                  IF (iflag(k,i).EQ.0) THEN
                    iflag(k,i) = 1
                    lun(k,i) = unit(j)
                    xmink(k,i) = xmin(j)
                    delk(k,i) = delx(j)
                    ireck(k,i) = irec(j)+npk*ntk*ngk*(ipo-1)
C                    print*,'k,i,ireck(k,i)',k,i,ireck(k,i)
                  ENDIF
               ENDIF
              ENDDO
C              print*,k,vwave(k),ireck(k,i),imatch
              IF(imatch.eq.0)THEN
                  print*,'Warning - read_klist. Match not found'
                  print*,'For wavelength : ',k,vwave(k)
                  print*,'Snapping to nearest wavelength'

                  MAXDX1=1000.0
                  DO ipo=1,npoint
                   dx = 100*ABS((vwave(k)-XCENK(ipo))/XCENK(ipo))
                   IF(dx.lt.MAXDX1)THEN
                    iflag(k,i) = 1
                    lun(k,i) = unit(j)
                    xmink(k,i) = xmin(j)
                    delk(k,i) = delx(j)
                    ireck(k,i) = irec(j)+npk*ntk*ngk*(ipo-1)
                    ipo1=ipo
                   ENDIF
                   MAXDX1 = dx
                  ENDDO
C                  print*,'read_klist',k,vwave(k),ipo1,xcenk(ipo1)

              ENDIF

             ENDDO
            ENDIF
          ENDIF
        ENDDO
30    ENDDO

C-----------------------------------------------------------------------
C
C	Finally, step through and and check that each wavenumber has
C	assigned variables for each gas. If not, write a warning and
C	flag the wavenumber/gas.
C
C-----------------------------------------------------------------------

      DO i=1,ngas
        DO j=1,nwave
          IF (iflag(j,i).EQ.0) THEN
C            WRITE(*,*)' READ_KLIST.f :: Error: Match not found in'
C            WRITE(*,*)' k-table for gas ',idgas(i),isogas(i)
C            WRITE(*,*)' at', vwave(j)
            lun(j,i) = -99
          ENDIF
        ENDDO
      ENDDO
				
C-----------------------------------------------------------------------
C
C	Return and end
C
C-----------------------------------------------------------------------

1020  FORMAT (A100)
1030  FORMAT (' Ktable filenames: ', A)
1035  FORMAT (' Calling read_khead. reading: ', A)

      RETURN

      END
