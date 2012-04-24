      subroutine openflux(fintname,ipop,LUNIS,nlayer,nmu,
     1 nwave,ng,nf,ioff,vwave,basep,baseh)
C     ********************************************************
C     Routine to open a binary internal radiation file for input
C     or output. If ipop=0, then file is opened for reading.
C     If ipop=1, then file is opened for writing
C     ********************************************************


      implicit none
      INCLUDE '../includes/arrdef.f'

      character*100 fintname
      integer nlayer,nmu,nwave,ng,nf,ioff,ipop,LUNIS,IRECL,I
      real vwave(maxbin),basep(maxlay),baseh(maxlay)

      irecl=4

      if(ipop.eq.0)then
       OPEN(LUNIS,FILE=fintname,STATUS='OLD',
     1    ACCESS='DIRECT',RECL=IRECL)

       READ(LUNIS,REC=1)nlayer
       READ(LUNIS,REC=2)nmu
       READ(LUNIS,REC=3)nwave
       READ(LUNIS,REC=4)ng
       READ(LUNIS,REC=5)nf

       DO I=1,nwave
        READ(LUNIS,REC=5+I)vwave(i)
       ENDDO

       DO I=1,nlayer
        READ(LUNIS,REC=5+NWAVE+I)basep(i)
       ENDDO

       DO I=1,nlayer
        READ(LUNIS,REC=5+NWAVE+NLAYER+I)baseh(i)
       ENDDO

       IOFF = 5+NWAVE+2*NLAYER+1

      ELSE

       OPEN(LUNIS,FILE=fintname,STATUS='UNKNOWN',
     1    ACCESS='DIRECT',RECL=IRECL)

       WRITE(LUNIS,REC=1)nlayer
       WRITE(LUNIS,REC=2)nmu
       WRITE(LUNIS,REC=3)nwave
       WRITE(LUNIS,REC=4)ng
       WRITE(LUNIS,REC=5)nf

       DO I=1,nwave
        WRITE(LUNIS,REC=5+I)vwave(i)
       ENDDO

       DO I=1,nlayer
        WRITE(LUNIS,REC=5+NWAVE+I)basep(i)
       ENDDO

       DO I=1,nlayer
        WRITE(LUNIS,REC=5+NWAVE+NLAYER+I)baseh(i)
       ENDDO

      ENDIF

      RETURN

      END      

