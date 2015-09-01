      REAL FUNCTION GETMASS(IDGAS,ISOGAS)
C     *************************************************************
C     Look up molecular mass of gas
C
C     Pat Irwin
C     16/2/12
C     *************************************************************

      IMPLICIT NONE
      INTEGER IDGAS,ISOGAS,J
      REAL XMASS

C------------------------------------------------------------------------------
C     Include line data base variables
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/dbcom.f'
C------------------------------------------------------------------------------

      IF(ISOGAS.EQ.0)THEN
C       Have to compute mass for normal combination if using all isotopes
        XMASS = 0.0
        DO 255 J=1,DBNISO(IDGAS)
          XMASS = XMASS + MASSNO(J,IDGAS)*RELABU(J,IDGAS)
255     CONTINUE

C       It is possible that the relative abundance for one of the isotopes is
C       wrong (e.g. set to 1 because unknown) so checking that the final mass is
C       within 20% of the mass of the main isotope
        IF(MASSNO(1,IDGAS).GT.1.E-32.AND.
     1  ABS((XMASS-MASSNO(1,IDGAS))/MASSNO(1,IDGAS)).GT.0.2)THEN
          XMASS = MASSNO(1,IDGAS)
c          WRITE(*,*)'getmass.f :: *WARNING* using main isotope'
c          WRITE(*,*)'mass-number.'
        ENDIF
      ELSE
        IF(ISOGAS.LE.DBNISO(IDGAS))THEN
          XMASS = MASSNO(ISOGAS,IDGAS)
        ELSE
          WRITE(*,*)'getmass.f :: *ERROR* Model doesn_t include'
          WRITE(*,*)'isotope ',ISOGAS,' for gas ',IDGAS,'.'
          WRITE(*,*)' '
          WRITE(*,*)'Stopping program.'
          STOP
        ENDIF
      ENDIF
 
      GETMASS=XMASS

      RETURN

      END
