      SUBROUTINE GWRITEPATL(RUNNAME,ISCAT,NCONV,VCONV,FWHM,LAYHT,NLAYER,
     1 LAYTYP,LAYINT,FLAGH2P)
C     $Id:
C     *******************************************************************
C     Subroutine to write out the .pat file needed for a CIRSradg run.
C
C     Input variables:
C       RUNNAME         character*100    Run root name
C       NCONV           integer         Number of convolution wavelengths
C       VCONV(MCONV)   	real 	        Convolution wavelengths
C	FWHM		real		FWHM of convolved spectrum
C       LAYHT           real            Base height of bottom layer 
C					 (if nadir)
C       NLAYER          integer         Number of layers
C       LAYTYP          integer         Type of layering
C       LAYINT          integer         Type of layer integration
C	FLAGH2P		integer		Equals 1 if para-H2 is variable
C
C     Pat Irwin	29/7/96		Original
C     Pat Irwin 17/10/03	Tidied for Nemesis
C
C     *******************************************************************

      IMPLICIT NONE
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      CHARACTER*100 RUNNAME
      REAL LAYHT,VCONV(MCONV),LAYANG,ILAYER,DELV,FWHM
      INTEGER LAYTYP,LAYINT,NLAYER,NCONV,LAYBOT,FLAGH2P,ISCAT
      CHARACTER*80 TEXT

      CALL FILE(RUNNAME,RUNNAME,'pat')

      OPEN(31,FILE=RUNNAME,STATUS='UNKNOWN')

      WRITE(31,1)' '
      TEXT='interval'
      WRITE(31,1)TEXT
1     FORMAT(A)
      IF(NCONV.GT.1)THEN
       DELV = (VCONV(NCONV)-VCONV(1))/FLOAT(NCONV-1)
      ELSE
       DELV=FWHM
      ENDIF 
      WRITE(31,*)VCONV(1),VCONV(NCONV),DELV,FWHM
      WRITE(31,1)'   24  0  0'
      WRITE(31,1)' '
 
      TEXT(1:10)='spec data '
      CALL FILE(RUNNAME,RUNNAME,'kls')
      TEXT(11:60)=RUNNAME(1:50)
      WRITE(31,1)TEXT
      WRITE(31,1)' '


      CALL FILE(RUNNAME,RUNNAME,'prf')
      TEXT=' '
      TEXT(1:6)='model '
      TEXT(7:60)=RUNNAME(1:54)
      WRITE(31,1)TEXT
      WRITE(31,*)' '

      TEXT = 'dust model aerosol.prf'
      WRITE(31,1)TEXT
      WRITE(31,*)' '


      CALL FILE(RUNNAME,RUNNAME,'xsc')
      TEXT=' '
      TEXT(1:13)='dust spectra '
      TEXT(14:60)=RUNNAME(1:47)
      WRITE(31,1)TEXT
      WRITE(31,*)' '

      IF(FLAGH2P.EQ.1)THEN
       TEXT = 'fparah2 model parah2.prf'
       WRITE(31,1)TEXT
       WRITE(31,*)' '
      ENDIF
          
2     format(a,i3)
3     format(a,f7.2)

      LAYANG = 90.0

      WRITE(31,1)'layer'
      WRITE(31,2)'nlay ',nlayer
      WRITE(31,3)'layht ',layht
      WRITE(31,3)'layang ',layang
      WRITE(31,2)'layint ',layint
      WRITE(31,2)'laytyp ',laytyp
      WRITE(31,1)' '

      DO 101 ILAYER=1,NLAYER

       WRITE(31,1)'atm'
       TEXT='limb'
       WRITE(31,5)TEXT,ILAYER
5      FORMAT(A6,I4)
       IF(ISCAT.EQ.0)THEN
        WRITE(31,1)'therm'
        WRITE(31,1)'noscatter'
       ELSE
        WRITE(31,1)'notherm'
        WRITE(31,1)'scatter'
       ENDIF
       WRITE(31,1)'nowf'
       WRITE(31,1)'nocg'
       WRITE(31,1)'noabsorb'
       WRITE(31,1)'binbb'
       WRITE(31,1)'nobroad'
       WRITE(31,1)' '

101   CONTINUE

      WRITE(31,1)'clrlay'
      WRITE(31,1)'nocombine'

      CLOSE(31)

      RETURN
      END



