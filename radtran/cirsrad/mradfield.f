      subroutine mradfield(imod) 
C     **************************************************************** 
C     Subroutine for modifying internal radiation field as per JPL 
C     MCS assumptions 
C     ****************************************************************
      implicit none
      include '../includes/arrdef.f'
      integer irecl,nlayerF,nmuF,nwaveF,ngF,lunis,imod
      integer i,ioff,iwave,ig,irec,j,k,nfF,i1
      real vwaveF(maxbin),basehf(maxlay),basepf(maxlay)
      real umif(1000,maxmu,maxscatlay,20,11)
      real uplf(1000,maxmu,maxscatlay,20,11),xtmp


      irecl=4
      lunis=61

      OPEN(LUNIS,FILE='internal.fld',STATUS='OLD',
     1    ACCESS='DIRECT',RECL=IRECL)

      READ(LUNIS,REC=1)nlayerF
      READ(LUNIS,REC=2)nmuF
      READ(LUNIS,REC=3)nwaveF
      READ(LUNIS,REC=4)ngF
      READ(LUNIS,REC=5)nfF


      DO I=1,nwaveF
          READ(LUNIS,REC=5+I)vwaveF(i)
      ENDDO

      DO I=1,nlayerF
          READ(LUNIS,REC=5+NWAVEF+I)basepF(i)
      ENDDO

      DO I=1,nlayerF
          READ(LUNIS,REC=5+NWAVEF+NLAYERF+I)basehF(i)
      ENDDO

      IOFF = 5+NWAVEF+2*NLAYERF+1

      do 70 iwave=1,nwavef
       do 65 ig=1,ngf

         irec = ioff+2*((iwave-1)*ngf + (ig-1))*nlayerf*nmuf*(nff+1)

         do j=1,nlayerf
          do k=1,nmuf
           do i1=1,nff+1
            read(LUNIS,REC=IREC)xtmp
            umif(iwave,k,j,ig,i1)=xtmp
            irec=irec+1
            read(LUNIS,REC=IREC)xtmp
            uplf(iwave,k,j,ig,i1)=xtmp
            irec=irec+1
           enddo
          enddo
         enddo

65     continue
70    continue
      close(LUNIS)


C     Now modify internal radiation field
      IF(IMOD.EQ.1)THEN
C      Set all upwelling radiation to that emerging straight upwards at
C      top of atmosphere and set all downwelling radiation to zero


       do 72 iwave=1,nwavef
        do 67 ig=1,ngf

         irec = ioff+2*((iwave-1)*ngf + (ig-1))*nlayerf*nmuf*(nff+1)

         do j=1,nlayerf
          do k=1,nmuf
           uplf(iwave,k,j,ig,1)=0.0
           umif(iwave,k,j,ig,1)=umif(iwave,nmuf,nlayerf,ig,1)
          enddo
         enddo

67      continue
72     continue

      ELSE IF(IMOD.EQ.2)THEN

C      Set all upwelling radiation to that emerging at 61.45 degrees at
C      top of atmosphere and set all downwelling radiation to zero


       do 73 iwave=1,nwavef
        do 68 ig=1,ngf

         irec = ioff+2*((iwave-1)*ngf + (ig-1))*nlayerf*nmuf*(nff+1)

         do j=1,nlayerf
          do k=1,nmuf
           uplf(iwave,k,j,ig,1)=0.0
           umif(iwave,k,j,ig,1)=umif(iwave,2,nlayerf,ig,1)
          enddo
         enddo

68      continue
73     continue

      ELSE
       print*,'modradfield: modification ID unrecognised',IMOD
       print*,'modradfield: doing nothing'
      ENDIF

C     Now write out modified radiation field

      OPEN(LUNIS,FILE='internal.fld',STATUS='UNKNOWN',
     1    ACCESS='DIRECT',RECL=IRECL)

      WRITE(LUNIS,REC=1)nlayerF
      WRITE(LUNIS,REC=2)nmuF
      WRITE(LUNIS,REC=3)nwaveF
      WRITE(LUNIS,REC=4)ngF
      WRITE(LUNIS,REC=5)nfF

      DO I=1,nwaveF
          WRITE(LUNIS,REC=5+I)vwaveF(i)
      ENDDO

      DO I=1,nlayerF
          WRITE(LUNIS,REC=5+NWAVEF+I)basepF(i)
      ENDDO

      DO I=1,nlayerF
          WRITE(LUNIS,REC=5+NWAVEF+NLAYERF+I)basehF(i)
      ENDDO

      IOFF = 5+NWAVEF+2*NLAYERF+1

      do 71 iwave=1,nwavef
       do 66 ig=1,ngf

         irec = ioff+2*((iwave-1)*ngf + (ig-1))*nlayerf*nmuf
         do j=1,nlayerf
          do k=1,nmuf
           do i1=1,nff+1
            xtmp = umif(iwave,k,j,ig,i1)
            write(LUNIS,REC=IREC)xtmp
            irec=irec+1
            xtmp=uplf(iwave,k,j,ig,i1)
            write(LUNIS,REC=IREC)xtmp
            irec=irec+1
           enddo
          enddo
         enddo

66     continue
71    continue

      close(LUNIS)

      return

      end
