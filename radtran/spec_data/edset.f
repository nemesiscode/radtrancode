      SUBROUTINE EDSET
C     $Id: edset.f,v 1.7 2011-06-23 09:09:02 irwin Exp $
C--------------------------------------------------------------
C
C_HIST:
C		NT	20.jan.2004 Added option to broaden PH3 for Saturn
C				h2 and he abundances.
C--------------------------------------------------------------
C     common variables used by all linedata routines
C     mixed data types are avoided
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
      LOGICAL WSBTAB,CO2TAB,H2OTAB,SCOFLD,WEDAD,BEZARD,PMIRR
      LOGICAL H2HePH3J,H2HePH3S,BERGH,BERGC,H2HeCH4
      COMMON /EDLOG/WSBTAB,CO2TAB,H2OTAB,SCOFLD,WEDAD,BEZARD,BERGH,
     $      BERGC,PMIRR,H2HePH3J,H2HePH3S,H2HeCH4
C--------------------------------------------------------------
C      character*2 QC(28)
C      REAL BH2(28),NH2(28),BHe(28),NHe(28)
C      INTEGER JN(28),KN(28)
C      COMMON /CH4EXOMOL/QC,BH2,NH2,BHe,NHe,JN,KN
C--------------------------------------------------------------



      LOGICAL ASKYN
      CHARACTER*256 BUFFER,DPNAME,IPNAME
1     FORMAT(A)
      WRITE(BUFFER,10)
10    FORMAT(' # generic EDLIN used')
      WRITE(3,1)BUFFER(1:DBRECL)

	write(*,*) '(answer no to all to use air broadened widths)'
      WSBTAB=ASKYN('set self broadened width to air broadened?')
      CO2TAB=ASKYN('multiply H2O air broadened width by 1.3?')
      H2OTAB=ASKYN('replace H2O air b. by delaye CO2-broadened data?')
      SCOFLD=ASKYN('generate linedata using PMR2 assumptions?')
      WEDAD=ASKYN('generate linedata using Wedads assumptions?')
      BEZARD=ASKYN('generate linedata using Brunos assumptions?')
      BERGH=ASKYN('generate linedata using Catherines assumptions?')
      BERGC=ASKYN('generate linedata using Bergh assump. (mod CH3D)?')
      PMIRR=ASKYN('generate linedata using all PMIRR assumptions?')
      H2HePH3J=ASKYN('H2-He broaden PH3 for Jupiter?')
      H2HePH3S=ASKYN('H2-He broaden PH3 for Saturn?')
      H2HeCH4=ASKYN('H2-He broaden CH4 from ExoMOL broad?')
      IF(WSBTAB)THEN
        WRITE(BUFFER,11)
11      FORMAT(' # setting self broadened widths to air broadened')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(CO2TAB)THEN
        WRITE(BUFFER,12)
12      FORMAT(' # multiplying H2O air broadened widths by 1.3')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(H2OTAB)THEN
        WRITE(BUFFER,13)
13      FORMAT(' # calculate H2O CO2 widths and CO2,H2O T dependencies')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(SCOFLD)THEN
        WRITE(BUFFER,14)
14      FORMAT(' # calculate H2O and CO2 widths/dependencies from JTS')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(WEDAD)THEN
        WRITE(BUFFER,15)
15      FORMAT(' # calculate H2O and CO2 widths/dependencies from 
     1 Wedad')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(BEZARD)THEN
        WRITE(BUFFER,16)
16     FORMAT(' # calculate line widths/dependencies from 
     1 Bruno Bezard Assumptions')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(BERGH)THEN
        WRITE(BUFFER,161)
161     FORMAT(' # calculate line widths/dependencies from 
     1 Catherine de Berghs Assumptions')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(BERGC)THEN
        WRITE(BUFFER,161)
171     FORMAT(' # calculate line widths/dependencies from 
     1 Catherine de Berghs Assump. (mod CH3D/CH4)')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(PMIRR)THEN
        WRITE(BUFFER,17)
17      FORMAT(' # H2O and CO2 widths/dependencies from 
     1 PMIRR assumptions')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(H2HePH3J)THEN
        WRITE(BUFFER,18)
18      FORMAT(' # H2-He broadened PH3 (Levy,A. etal 1993) for Jupiter')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(H2HePH3S)THEN
        WRITE(BUFFER,19)
19      FORMAT(' # H2-He broadened PH3 (Levy,A. etal 1993) for Saturn')
        WRITE(3,1)BUFFER(1:DBRECL)
      END IF
      IF(H2HeCH4)THEN
        WRITE(BUFFER,19)
20      FORMAT(' # H2-He broadened CH4 from ExoMOL')
        WRITE(3,1)BUFFER(1:DBRECL)
 
C        DPNAME='/network/group/aopp/planetary/PGJI001_IRWIN_LBLKTAB'
C        DPNAME=DPNAME//'/linedata/hitran20/EXOMOL-CH4'
C        IPNAME=DPNAME//'12C-1H4__H2.broad.txt'
C        OPEN(12,FILE=IPNAME,status='old')
C         DO I=1,11
C          READ(12,101)QC(I),BH2(I),NH2(I),JN(I),KN(I)
C         ENDDO
C         DO I=12,28
C          READ(12,102)QC(I),BH2(I),NH2(I),JN(I)
C         ENDDO
C        CLOSE(12)
C        IPNAME=DPNAME//'12C-1H4__He.broad.txt'
C        OPEN(12,FILE=IPNAME,status='old')
C         DO I=1,11
C          READ(12,101)QC(I),BHe(I),NHe(I),JN(I),KN(I)
C         ENDDO
C         DO I=12,28
C          READ(12,102)QC(I),BHe(I),NHe(I),JN(I)
C         ENDDO
C        CLOSE(12)
      END IF

C101   FORMAT(A2,F6.4,F5.3,I7,I2)
C102   FORMAT(A2,F6.4,F5.3,I7)

      RETURN
      END
