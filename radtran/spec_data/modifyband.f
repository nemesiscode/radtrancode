      program modifyband
      integer iread(10),IP(10),idgas(10),isogas(10),ij1(10),ij2(10)
      real twaven(3000,2),tpout(3000,10,7)
      real swaven(3000,2),spout(3000,10,7)
      logical old

      character*100 ipfile,opfile
      character*20 head
      character*10 buffer

      print*,'Enter name of band file : '
      read(5,1)ipfile

      call file(ipfile,ipfile,'ban')

      open(12,file=ipfile,status='old')
c     first skip header     
11    read(12,500)buffer
      if(buffer.ne.'**********')goto11

      read(12,401)vmin
      read(12,402)delv
      read(12,400)fwhm
      read(12,403)npoint
      read(12,404)ngas
      print*,'ngas = ',ngas
      do i=1,ngas
       read(12,405)idgas(i),isogas(i)
       print*,idgas(i),isogas(i)
      enddo

      read(12,406)head
      
      do 105 i=1,ngas
       READ(12,406)HEAD
       READ(12,406)HEAD
       IF(HEAD(1:3).EQ.'Ex.')THEN
        DO 111 J=1,NPOINT
         READ(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         READ(12,933)TPOUT(J,I,6),TPOUT(J,I,7)
111     CONTINUE
       ELSE
        DO 110 J=1,NPOINT
         READ(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         TPOUT(J,I,6)=0.0
         TPOUT(J,I,7)=0.0
110     CONTINUE
       ENDIF
105   continue

      close(12)


      vmax = vmin + (npoint-1)*delv
      
      call prompt('Enter output filename : ')
      read(5,1)opfile

      call file(opfile,opfile,'ban')
      

      open(12,file=opfile,status='new')
      write(12,401)vmin
      write(12,402)delv
      write(12,400)fwhm
      write(12,403)npoint
      write(12,404)ngas
      do i=1,ngas
       write(12,405)idgas(i),isogas(i)
      enddo

      write(12,406)'---------------------------------------------------------
     &---------------------'
      
      do 305 I=1,ngas
       WRITE(12,410)IDGAS(I),ISOGAS(I)
       WRITE(12,406)'Ex. V0       dV       Kv(T0)    delta/AD0
     &      y0           El       SFB'
        
       DO 311 J=1,NPOINT

         IF(TWAVEN(J,1).GE.4000.0.AND.TWAVEN(J,1).LE.11000.0)THEN
          TPOUT(J,I,5)=5.5
         ENDIF
         WRITE(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         WRITE(12,933)TPOUT(J,I,6),TPOUT(J,I,7)

311    CONTINUE
305   continue

      close(12)
     

500   FORMAT(1X,A10)
401   FORMAT(1X,'VMIN = ',F8.2)
402   FORMAT(1X,'DELV = ',F8.2)
400   FORMAT(1X,'FWHM = ',F8.2)
403   FORMAT(1X,'NPOINT = ',I5)
404   FORMAT(1X,'NGAS = ',I3)      
405   FORMAT(1X,'Gas : ',I3,I3)
406   FORMAT(1X,A)
410   FORMAT(1X,'Band data for gas : ',I3,I3)
932   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5)
933   FORMAT(1x,2(e12.5))
1     format(a)

      END
