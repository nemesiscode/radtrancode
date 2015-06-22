      SUBROUTINE RESERVEGAS
C     **************************************************************
C     Subroutine to read in relative abundances of gases from the
C     reference gasinfo.dat file, which is expected to reside in the
C     raddata/ directory.
C
C     Pat Irwin		Original	9/5/12
C
C     **************************************************************
      INCLUDE '../includes/dbcom.f'

      GASFIL='gasinforef.dat'
      CALL DATARCHIVE(GASFIL)
      DBLUN=2

      CALL RDGAS

      RETURN

      END
      
