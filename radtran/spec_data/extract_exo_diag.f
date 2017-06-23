      PROGRAM EXTRACT_EXO_DIAG
C     $Id:
C---------------------------------------------------------------------------
C
C_TITLE: EXTRACT_EXO_DIAG - extract band diagnostic information from ExoMOL
C
C_ARGS:  None.
C
C_KEYS:  PROG, VMS, CAL .
C
C_DESCR: Prompts for line data key and wavelength range.
C
C_FILES: unit DBLUN  line data and gas files
C
C_CALLS:  LINESS        reads line data
C         RDKEY         reads in details from the key file
C         PROMPT        user prompt utility
C         REMSP         removes leading spaces from text string
C         RDGAS         reads in gas information
C         RDISO         reads in isotope information
C        			
C_BUGS:  
C
C_HIST:  23jun17 PGJI ORIGINAL VERSION adapted from li_lines.f
C
C_END:
C
C--------------------------------------------------------------------------
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      INTEGER LINLIM
      PARAMETER (LINLIM=200000)
      REAL VLIN(LINLIM),ALIN(LINLIM),ELIN(LINLIM),
     1 SBLIN(LINLIM),TDW(LINLIM),TDWS(LINLIM),PSHIFT(LINLIM),
     2 DOUBV(LINLIM)
      DOUBLE PRECISION SLIN(LINLIM),sums,sume,sumaw,sumsw
      CHARACTER*15 LLQ(LINLIM)
      INTEGER IDLIN(LINLIM),IDGAS(MAXDGAS),ISOGAS(MAXDGAS)
      INTEGER NGAS,IGAS,MAXLIN,IFIRST,NLIN,ILAST,I,DISEXP
      INTEGER NEXTRA
      CHARACTER*1 ANS

      REAL VMIN,VMAX,DELV
      CHARACTER*12 MANTIS
C
      CALL PROMPT('data base key?')
      READ(*,102)KEYFIL
      CALL REMSP(KEYFIL)
102   FORMAT(A)
      CALL RDKEY(2)
      CALL RDGAS
      CALL RDISO
      NGAS=1
      IDGAS(1)=11
      ISOGAS(1)=0

      print*,'ngas,id,iso : ',NGAS,IDGAS(1),ISOGAS(1)
      CALL PROMPT('Enter wavelength minumum and maximum : ')
      READ*,X1,X2
 
      CALL PROMPT('Enter step and FWHM : ')
      READ*,XDELV,XFWHM

      NX = 1 + INT((X2-X1)/XDELV)



      OPEN(12,file='extract.out',status='unknown')
      WRITE(12,*)X1,X2,XDELV,XFWHM
      
      DO 111 I=1,NX
       x = X1+(I-1)*XDELV
       xa = x-0.5*XFWHM
       xb = xa+XFWHM
       vmin = 1e4/xb
       vmax = 1e4/xa
       print*,I,vmin,vmax
       MAXLIN=LINLIM
       NEXTRA=LINLIM
       IFIRST=1
       DELV=VMAX-VMIN

       CALL LINESS(VMIN,DELV,MAXLIN,VLIN,SLIN,ALIN,ELIN,IDLIN,SBLIN,
     1PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,IFIRST,ILAST,NGAS,IDGAS,ISOGAS)

       WRITE(*,110)NLIN
110    FORMAT(' found ',I8,' lines')

       sums=0.
       sume=0.
       sumaw=0.
       sumsw=0.
       do j=1,nlin
C        print*,slin(j),elin(j),alin(j),sblin(j)
        sums=sums+slin(j)
        sume=sume+dble(elin(j)*slin(j))
        sumaw=sumaw+dble(alin(j)*slin(j))
        sumsw=sumsw+dble((alin(j)-sblin(j))*slin(j))
       enddo
       if(nlin.gt.0)then
        sume=sume/sums
        sumaw=sumaw/sums
        sumsw=sumsw/sums

        print*,X,vmin,vmax,nlin,sums,sums/float(nlin),sume,sumaw,sumsw
        WRITE(12,*)X,vmin,vmax,nlin,sums,sums/float(nlin),sume,sumaw,
     1   sumsw
       else
        print*,X,vmin,vmax,nlin,sums,sums,sume,sumaw,sumsw
        WRITE(12,*)X,vmin,vmax,nlin,sums,sums,sume,sumaw,sumsw 
       endif
111   CONTINUE
      CLOSE(12)

      END
