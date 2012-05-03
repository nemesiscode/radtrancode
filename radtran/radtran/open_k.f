	SUBROUTINE OPEN_K(LUN1, XMIN, nx, delx, fwhm, IDGAS,
     1		ISOGAS, NG, np, nt, press1, temp1, del_g, g_ord, irec)

C     ****************************************************************
C     Subroutine to open a .kta k-coefficient lookup file for output
C
C     Input variables
C	LUN1		INTEGER	logical unit number
C	XMIN		REAL	Minimum wavelength/wavenumber
C       NX              INTEGER Number of bins
C       DELX		REAL    Seperation of bin centers
C	IDGAS		INTEGER Radtran gas identifier
C	ISOGAS		INTEGER Radtran isotope identifier
C	NG		INTEGER	Number of g-ordinates
C	NP		INTEGER	Number of pressure points
C	NT		INTEGER Number of temperature points
C	PRESS1(20)	REAL	Pressure grid
C	TEMP1(20)	REAL	Temperature grid
C	DEL_G(NG)	REAL	Gaussian g-space weights
C	G_ORD(NG)	REAL	Gaussian g-space ordinates
C
C     Output variables
C	IREC		INTEGER	Index of first writeable record
C
C     Pat Irwin		23/2/96
C
C     ****************************************************************

	implicit none
	integer		NP, NT, NG, IRECL, LUN1, IDGAS, ISOGAS,
     1			IREC0, nx, J, ISYS, IREC
      	real 		PRESS1(20),TEMP1(20),G_ORD(ng),DEL_G(ng),
     1			delx, fwhm, xmin
	character*100	KTAFIL, OPFILE 

      write (*,'(''Enter output filename: '',$)')
      READ(5,1)OPFILE
      CALL FILE(OPFILE,KTAFIL,'kta')
1     FORMAT(A)

      IRECL=ISYS()			! New tables use 1 words/record
      OPEN(UNIT=LUN1,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1RECL=IRECL)
      IREC0=11 + 2*NG + 2 + NP+NT + 2
      WRITE(LUN1,REC=1)IREC0
      WRITE(LUN1,REC=2)nx
      WRITE(LUN1,REC=3)xmin
      WRITE(LUN1,REC=4)delx
      WRITE(LUN1,REC=5)fwhm
      WRITE(LUN1,REC=6)NP
      WRITE(LUN1,REC=7)NT
      WRITE(LUN1,REC=8)NG
      WRITE(LUN1,REC=9)idgas
      WRITE(LUN1,REC=10)isogas
      IREC=11
      DO 299 J=1,NG
       WRITE(LUN1,REC=IREC)G_ORD(J)
       IREC=IREC+1
299   CONTINUE
      DO 399 J=1,NG
       WRITE(LUN1,REC=IREC)DEL_G(J)
       IREC=IREC+1
399   CONTINUE
      IREC=11 + 2*NG + 2
      DO 301 J=1,NP
       WRITE(LUN1,REC=IREC) PRESS1(J)
       IREC=IREC+1
301   CONTINUE
      DO 302 J=1,NT
       WRITE(LUN1,REC=IREC)TEMP1(J)
       IREC=IREC+1
302   CONTINUE


      IREC=IREC0
      RETURN
      END

