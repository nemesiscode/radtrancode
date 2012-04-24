      SUBROUTINE WTEXT(TEXT)
C     $Id: wtext.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ***********************************************************************
C     Simply writes text to the screen
C     Seperate program so that carriage control can be varied between systems
C     VMS/UNIX version
C
C	1/1/90	SBC	Original Version
C	3/10/94	PGJI	Added $Id: wtext.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C
C     ***********************************************************************
      CHARACTER TEXT*(*) 
      WRITE(*,10)TEXT
10    FORMAT(A)
      RETURN
      END
