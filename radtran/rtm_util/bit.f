      LOGICAL FUNCTION BIT(I,FLAG)
C     $Id: bit.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ********************************************************************
C     Returns .true. if bit I of FLAG is 1
C     Lowest significant bit is bit zero
C
C     1/1/90    Original Version:       SBC
C     3/10/94   Updated Header          PGJI
C
C     ********************************************************************

      INTEGER I,FLAG
      IF( (FLAG/(2**I)) .EQ. (2*(FLAG/(2**(I+1)))) )THEN
        BIT=.FALSE.
       ELSE
        BIT=.TRUE.
        END IF
      RETURN
      END
