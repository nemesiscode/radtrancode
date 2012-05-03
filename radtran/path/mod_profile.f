      PROGRAM MOD_PROFILE
C     $Id:
C     *********************************************************************
C     Modifies specific density values in dust profile by user-specified
C     fixed value.
C
C     Pat Irwin    6/2/02
C
C     *********************************************************************
      IMPLICIT NONE
      include '../includes/arrdef.f'
      CHARACTER*100 IPFILE,OPFILE
      CHARACTER*100 BUFFER
      INTEGER NPRO,NCONT,I,J
      REAL CONT(MAXCON),XCORR(MAXCON),H
C     *******************************************


      CALL PROMPT('Enter name of input dust profile?')
      READ(*,10)IPFILE
10    FORMAT(A)
      CALL FILE(IPFILE,IPFILE,'prf')
      CALL PROMPT('Enter name of output modified dust profile?')
      READ(*,10)OPFILE
      CALL FILE(OPFILE,OPFILE,'prf')


      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      OPEN(UNIT=2,FILE=OPFILE,STATUS='UNKNOWN')

      READ(1,10)BUFFER

      BUFFER='# '
      DO I=1,LEN(OPFILE)
       BUFFER(I+2:I+2)=OPFILE(I:I)
      END DO

      WRITE(2,10)BUFFER
      READ(1,*)NPRO, NCONT
      WRITE(2,*)NPRO, NCONT

      PRINT*,'There are ',NCONT,' different particle types'
      CALL PROMPT('Enter OD correction factors : ')
      READ*,(XCORR(I),I=1,NCONT)

      DO 31 I=1,NPRO
        READ(1,*) H,(CONT(J),J=1,NCONT)
        DO J=1,NCONT
         CONT(J)=CONT(J)*XCORR(J)
        ENDDO
        WRITE(2,*) H,(CONT(J),J=1,NCONT)
31    CONTINUE
      CLOSE(1)
      CLOSE(2)
      END


