C***********************************************************************
C_TITL:	READ_KLBLLIST
C
C_DESC:	Writes variables to common block INTERPKLBL for subsequent use in 
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
C_FILE:	unit=11		klist (<runname>.lls)
C
C_CALL:	FILE		Forces file extension.
C	READ_KLBLHEAD	Opens .lta file, reads in header and returns
C			the beginning record number for the user-selected
C			spectral range as well as temp and pressure ranges
C			used.
C	FINDLOC		For a given monotonic array and 2 values, finds
C			the indices of the array values contained by the
C			passed variables.
C
C_HIST:	???????	???	Orginal version.
C***********************************************************************

      SUBROUTINE read_klbllist (klist,ngas,idgas,isogas,nwave,vwave)

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
C Defines the maximum values for a series of variables (layers, bins,
C paths, etc.)

      INTEGER maxkfil,imatch,ipo1,jpo
      real MAXDX,MAXDX1,ta1,ta2
      PARAMETER (maxkfil=100,MAXDX=1e-4)

      INTEGER ngas,nwave,nkl,npk,ntk,i,j,i1,i2,k,npointk
C NPK: Number of pressure points in k-tables.
C NTK: Number of temperature points in k-tables.
      INTEGER idgas(ngas),isogas(ngas),lun(maxgas)
      INTEGER unit(maxkfil),npt(maxkfil),idgask(maxkfil),
     1 isogask(maxkfil),np(maxkfil),nt(maxkfil),irec(maxkfil)
      INTEGER irec0,ipo

      REAL pmin,pmax,tmin,tmax,delvk
      REAL vwave(nwave)
      REAL xmin(maxkfil),xmax(maxkfil),delx(maxkfil)
      REAL pk(maxk),tk(maxk),t2k(maxk,maxk)
C PK: Pressure values in k-tables.
C TK: Temperature values in k-tables.
      REAL xmink
C XMINK: Beginning wavenumber in table.
C DELK: Wavenumber step of evaluation points.
      REAL kout(maxlay,maxgas),dkoutdt(maxlay,maxgas)
      CHARACTER*100 klist
      CHARACTER*200 ktafil(maxkfil),null

      COMMON /interpklbl/ lun,irec0,xmink,delvk,npointk,pk,npk,tk,
     1 t2k,ntk,kout,dkoutdt

C********************************* CODE ********************************

C-----------------------------------------------------------------------
C
C	Read in names of files.
C
C-----------------------------------------------------------------------

      nkl = 0
      null = ' '
      print*,'klbllist file : ',klist
      OPEN (UNIT=11,FILE=klist,STATUS='old')

10    nkl = nkl + 1
      READ(11,1020,END=20)ktafil(nkl)
      IF(ktafil(nkl).EQ.null)GOTO 20

      CALL file(ktafil(nkl),ktafil(nkl),'lta')
      WRITE(*,1030)ktafil(nkl)
      GOTO 10
	
20    nkl = nkl - 1
      CLOSE(UNIT=11)

C-----------------------------------------------------------------------
C
C	Open appropriate LBL lookup files and compare. 
C
C-----------------------------------------------------------------------
      WRITE(*,*)' READ_KLBLLIST.f :: Number of l-files = ',nkl

      DO i=1,nkl
        unit(i) = 100 + i
        WRITE(*,1035)KTAFIL(i)
        CALL read_klblhead(ktafil(i),unit(i),npt(i),xmin(i),delx(i),
     1   idgask(i),isogask(i),pk,tk,t2k,np(i),nt(i),irec(i))
        xmax(i) = xmin(i) + (npt(i) - 1) * delx(i)

        IF (i.EQ.1) THEN
          pmin = pk(1)
          pmax = pk(np(i))
          if(nt(i).gt.0)then
           tmin = tk(1)
           tmax = tk(nt(i))
          else
           tmin = t2k(1,1)
           tmax = t2k(np(i),abs(nt(i)))
          endif
          npointk = npt(i)
          xmink = xmin(i)
          delvk = delx(i)
          irec0 = irec(i)
          if(npointk.gt.mpoint)then
           print*,'Error in READ_KLBLLIST: npoint > mpoint'
           print*,npointk,mpoint
           print*,'Use smaller k-tables or modify MPOINT in arrdef.f'
           print*,'and recompile everything'
           stop
          endif
          npk = np(i)
          ntk = nt(i)
          delvk = delx(i)
        ELSE
          IF((np(i).NE.npk).OR.(nt(i).NE.ntk).OR.
     &		(delx(i).NE.delvk).OR.(npt(i).NE.npointk).OR.
     &           (xmin(i).NE.xmink).OR.(irec0.NE.irec(i))) THEN
            WRITE(*,*)' READ_KLBLLIST.f :: Error: Problems reading'
            WRITE(*,*)' lbl-tables array sizes. Stopping program.'
            WRITE(*,*)' '
            WRITE(*,*)' ktafil(i) ',ktafil(i)
            WRITE(*,*)' npk = ',npk,' np(i) = ',np(i)
            WRITE(*,*)' ntk = ',ntk,' nt(i) = ',nt(i)
            WRITE(*,*)' delvk = ',delvk,' delx(i) = ',delx(i)
            WRITE(*,*)' npointk = ',npointk,' npt(i) = ',delx(i)
            WRITE(*,*)' xmink = ',xmink,' xmin(i) = ',xmin(i)
            WRITE(*,*)' irec0 = ',irec0,' irec(i) = ',irec(i)
            STOP
          ENDIF

          if(nt(i).gt.0)then
           ta1 = tk(1)
           ta2 = tk(nt(i))
          else
           ta1 = t2k(1,1)
           ta2 = t2k(np(i),abs(nt(i)))
          endif

          IF ((abs(pk(1)-pmin).GT.MAXDX).OR.
     1    (abs(pk(npk)-pmax).GT.MAXDX).OR.
     2    (abs(ta1-tmin).GT.MAXDX).or.
     3    (abs(ta2-tmax).GT.MAXDX)) THEN
            WRITE(*,*)' READ_KLBLLIST.f :: Error: Problems reading'
            WRITE(*,*)' lbl-tables PT Grid. Stopping program.'
            WRITE(*,*)' '
            WRITE(*,*)' ktafil(I) ',ktafil(I)
            WRITE(*,*)' pmin = ',pmin,' pk(1) = ',pk(1)
            WRITE(*,*)' pmax = ',pmax,' pk(npk) = ',pk(npk)
            WRITE(*,*)' tmin = ',tmin,' ta1 = ',ta1
            WRITE(*,*)' tmax = ',tmax,' ta2 = ',ta2
            STOP
          ENDIF

        ENDIF
      ENDDO
	
C-----------------------------------------------------------------------
C
C	Check some variable sizes.
C
C-----------------------------------------------------------------------

      IF ((npk.GT.maxk).OR.(abs(ntk).GT.maxk)) THEN
        WRITE(*,*)' READ_KLBLLIST.f :: Error: Too many P/T points in'
        WRITE(*,*)' k-tables. Stopping program.' 
        WRITE(*,*)' '
        WRITE(*,*)' maxk = ',maxk
        WRITE(*,*)' np = ',np,' nt = ',nt
        STOP
      ENDIF


C     See which gases are included and assign unit numbers
      DO i=1,ngas
        lun(i)=-99
        DO j=1,nkl
          IF ((idgas(i).EQ.idgask(j)).AND.
     1    (isogas(i).EQ.isogask(j))) THEN
              lun(i)=unit(j)
          ENDIF
        ENDDO
      ENDDO
				
C-----------------------------------------------------------------------
C
C	Return and end
C
C-----------------------------------------------------------------------

1020  FORMAT (A200)
1030  FORMAT (' LBLtable filenames: ', A)
1035  FORMAT (' Calling read_klblhead. reading: ', A)

      RETURN

      END
