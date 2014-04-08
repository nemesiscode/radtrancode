      PROGRAM PARAH2_PROFILE
C     $Id%
C***********************************************************************
C_TITL:	PARAH2_PROFILE.f
C
C_DESC:	Generates a para-H2 profile file for use by RADTRANS where
C       where the para-H2 fraction is specified as a function of height
C	If the fraction is set negative, then the para-H2 fraction is
C       assumed to be in thermal equilibrium
C
C_ARGS: See definitions below.
C
C_FILE:	unit=1		<ipfile>.prf, and then again for <outfile>.prf
C
C_CALL:	prompt		Prompts the user for input.
C	file		Forces file extension.
C
C_HIST:	23/10/03	PGJI	ORIGINAL VERSION
C***********************************************************************

      IMPLICIT NONE

      INTEGER i,j,k,npro,nvmr,n,amform,iplanet
C IPLANET: Planet identification code (3= Earth, 5= Jupiter, etc.).

      INTEGER gasid(15),gasiso(15)
C GASID: Radtran gas identification code.
C GASISO: Radtran gas-isotope identification code.

      REAL molwt,radius,fpara,calcpara,xpara
C MOLWT: Molecular weight [kmol/kg].
C RADIUS: Planetary radius [km] at given latitude.

      REAL h(1000),p(1000),t(1000),vmr(1000,15)
C H: Height above reference surface [km].
C P: Pressure profile [atm].
C T: Temperature [K].


      REAL d,latitude
C LATITUDE: Latitude [degrees].

      CHARACTER*100 ipfile,outfile
C IPFILE: Input profile name.
C OUTFILE: Output profile name.
      CHARACTER*100 buffer

C********************************* CODE ********************************

C-----------------------------------------------------------------------
C 
C	Read in and modify the existing T/P/vmr profile to incorporate
C	the condensed vmr data.
C
C	The nominal format (AMFORM=0) for .ref/.prf files are:
C	AMFORM
C	IPLANET,LAT,NPRO,NVMR,MOLWT
C	GASID(1),GASISO(1)
C	:        :
C	:        :
C	GASID(NVMR),GASISO(NVMR)
C	*** COLUMN HEADER ***
C	H(1)    P(1)    T(1)    VMR(1,1)    VMR(1,2)    VMR(1,3)
C	:       :       :       :           :           :
C	:       :       :       :           :           :
C	H(NPRO) P(NPRO) T(NPRO) VMR(NPRO,1) VMR(NPRO,2) VMR(NPRO,3)
C	*** CONTINUED COLUMN HEADER (IF NECESSARY) ***
C	VMR(1,4)  VMR(1,5)      ...     VMR(1,NVMR)
C	:         :                     :
C	:         :                     :
C	VMR(NP,4) VMR(NP,5)     ...     VMR(NP,NVMR)
C
C	i.e. profiles are stored in blocks of 6 columns each with a
C	descriptive header. AMFORM=1 assumes that vmrs add up to 1.0
C       and so molecular weight can be calculated at each level.
C 
C-----------------------------------------------------------------------  

      WRITE(*,*)'Enter name of temp/press profile file: '
      READ(*,10)ipfile
10    FORMAT(A)
      CALL FILE(ipfile,ipfile,'prf')
      OPEN(UNIT=1,FILE=ipfile,STATUS='OLD')

C     Skip the header
54    READ(1,10)buffer
      IF(buffer(1:1).EQ.'#')GOTO 54
      READ(buffer,*)amform
      IF(AMFORM.EQ.1)THEN
       READ(1,*)iplanet,latitude,npro,nvmr
      ELSE
       READ(1,*)iplanet,latitude,npro,nvmr,molwt
      ENDIF
C Read in the gas-identificaiton information
      DO 20 i=1,nvmr
        READ(1,*)gasid(i),gasiso(i)
20    CONTINUE

C Skip header
      READ(1,*)
      DO 30 i=1,npro
          READ(1,*)h(i),p(i),t(i),(vmr(i,j),j=1,nvmr)
30    CONTINUE

      CLOSE(UNIT=1)


C***********************************************************************

      WRITE(*,*)' Enter output name of para-H2 file: '
      READ(*,10)outfile

      print*,'Enter para-H2 fraction'
      print*,'0 = calculate eqm fraction for given temperature'
      CALL PROMPT('-1 = always set to eqm) : ') 

      READ*,FPARA

C Write the output, without header to OUTFILE.prf
      CALL FILE(outfile,outfile,'prf')
      OPEN(UNIT=1,FILE=outfile,STATUS='UNKNOWN')
      WRITE(1,*)npro
      DO 31 I=1,npro
        xpara = fpara
        if(xpara.eq.0.0)then
         xpara = calcpara(dble(t(i)))
         print*,t(i),xpara
        endif
        WRITE(1,*)h(i),xpara
31    CONTINUE

      CLOSE(UNIT=1)

      END
