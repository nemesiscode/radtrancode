      program mergeband
      integer idgas1(20),isogas1(20),idgas2(20),isogas2(20)
      integer iread(20),IP(20),idgas(20),isogas(20),ij1(20),ij2(20)
      real twaven(3000,2),tpout(3000,20,7)
      real swaven(3000,2),spout(3000,20,7)
      logical old

      character*100 ipfile1,ipfile2,opfile
      character*20 head
      character*10 buffer

      print*,'Enter names of two band files : '
      read(5,1)ipfile1
      read(5,1)ipfile2

      call file(ipfile1,ipfile1,'ban')

      open(12,file=ipfile1,status='old')
c     first skip header     
11    read(12,500)buffer
      if(buffer.ne.'**********')goto11

      read(12,401)vmin1
      read(12,402)delv1
      read(12,400)fwhm1
      read(12,403)npoint1
      read(12,404)ngas1
      print*,'ngas1 = ',ngas1
      do i=1,ngas1
       read(12,405)idgas1(i),isogas1(i)
       print*,idgas1(i),isogas1(i)
      enddo

      read(12,406)head
      
      do 105 i=1,ngas1
       READ(12,406)HEAD
       READ(12,406)HEAD
       IF(HEAD(1:3).EQ.'Ex.')THEN
        DO 111 J=1,NPOINT1
         READ(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         READ(12,933)TPOUT(J,I,6),TPOUT(J,I,7)
111     CONTINUE
       ELSE
        DO 110 J=1,NPOINT1
         READ(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         TPOUT(J,I,6)=0.0
         TPOUT(J,I,7)=0.0
110     CONTINUE
       ENDIF
105   continue

      close(12)



      call file(ipfile2,ipfile2,'ban')

      open(12,file=ipfile2,status='old')
c     first skip header     
21    read(12,500)buffer
      if(buffer.ne.'**********')goto21

      read(12,401)vmin2
      read(12,402)delv2
      read(12,400)fwhm2
      read(12,403)npoint2
      read(12,404)ngas2
      print*,'ngas2 = ',ngas2
      do i=1,ngas2
       read(12,405)idgas2(i),isogas2(i)
       print*,idgas2(i),isogas2(i)
      enddo

      read(12,406)head
      
      do 205 i=1,ngas2
       READ(12,406)HEAD
       READ(12,406)HEAD
       IF(HEAD(1:3).EQ.'Ex.')THEN
        DO 211 J=1,NPOINT2
         READ(12,932)SWAVEN(J,1),SWAVEN(J,2),SPOUT(J,I,1),SPOUT(J,I,2),
     1 SPOUT(J,I,3),SPOUT(J,I,4),SPOUT(J,I,5)
         READ(12,933)SPOUT(J,I,6),SPOUT(J,I,7)
211     CONTINUE
       ELSE
        DO 210 J=1,NPOINT2
         READ(12,932)SWAVEN(J,1),SWAVEN(J,2),SPOUT(J,I,1),SPOUT(J,I,2),
     1 SPOUT(J,I,3),SPOUT(J,I,4),SPOUT(J,I,5)
         TPOUT(J,I,6)=0.0
         TPOUT(J,I,7)=0.0
210     CONTINUE
       ENDIF
205   continue

      close(12)

      if(delv1.ne.delv2)then
       print*,'delv1 <> delv2',delv1,delv2
       stop
      endif

      if(fwhm1.ne.fwhm2)then
       print*,'fwhm1 <> fwhm2',fwhm1,fwhm2
       stop
      endif

      vmax1 = vmin1 + (npoint1-1)*delv1
      vmax2 = vmin2 + (npoint2-1)*delv2

      vmin = min(vmin1,vmin2)
      vmax = max(vmax1,vmax2)
      delv = delv1
      fwhm = fwhm1

      do i=1,20
       ij1(i)=-1
       ij2(i)=-1
      end do

      do i=1,ngas1
       ico=i
       idgas(i)=idgas1(i)
       isogas(i)=isogas1(i)
       iread(i)=1
       ij1(i)=i
      enddo
      
      do i=1,ngas2
       old = .false.
       do j=1,ico
        if((idgas2(i).eq.idgas(j)).and.(isogas2(i).eq.isogas(j)))then
         old=.true.
         jold=j
        endif
       enddo
       if(old)then
        iread(jold)=3
        ij2(jold)=i
       else
        ico=ico+1
        idgas(ico) = idgas2(i)
        isogas(ico) = isogas2(i)
        iread(ico) = 2
        ij2(ico)=i
       endif
      enddo

      ngas = ico
      print*,'Ngas = ',ngas
      print*,'ID   ISO   IFILE   IJ1    IJ2'
      do i=1,ngas
       print*,idgas(i),isogas(i),iread(i),ij1(i),ij2(i)
       ip(i)=iread(i)
       if(iread(i).eq.3)then
        print*,'Two data sources for this gas. Enter priority 
     &file number : '
        read*,ip(i)
       endif
      end do

      npoint = 1 + int((vmax-vmin)/delv)

      print*,'vmin1,vmax1,npoint1 = ',vmin1,vmax1,npoint1
      print*,'vmin2,vmax2,npoint2 = ',vmin2,vmax2,npoint2
      print*,'vmin,vmax,npoint = ',vmin,vmax,npoint
      
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
      
      do 305 k=1,ngas
       WRITE(12,410)IDGAS(k),ISOGAS(k)
       WRITE(12,406)'Ex. V0       dV       Kv(T0)    delta/AD0
     &    y0      El       SFB'
        
       print*,k,idgas(i),isogas(k)
       DO 311 JJ=1,NPOINT
        VV = VMIN + (JJ-1)*DELV
        IS1 = -1
        IS2 = -1
        IF(VV.GE.VMIN1.AND.VV.LE.VMAX1)IS1 = 1
        IF(VV.GE.VMIN2.AND.VV.LE.VMAX2)IS2 = 1
        J1 = 1 + INT((VV - VMIN1)/DELV1)
        J2 = 1 + INT((VV - VMIN2)/DELV2)
        IF(J1.LT.0)J1=-1
        IF(J1.GT.NPOINT1)J1=-1
        IF(J2.LT.0)J2=-1
        IF(J2.GT.NPOINT2)J2=-1


        IF(IS1.EQ.1.AND.IS2.EQ.-1)IDS = 1
        IF(IS2.EQ.1.AND.IS1.EQ.-1)IDS = 2
        IF(IS1.EQ.1.AND.IS2.EQ.1)IDS = IP(K)

        IF(IDS.EQ.1)THEN
          I = IJ1(K)
          J = J1
          WRITE(12,932)VV,FWHM,TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
          WRITE(12,933)TPOUT(J,I,6),TPOUT(J,I,7)
        ELSE
          I = IJ2(K)
          J=J2
          WRITE(12,932)VV,FWHM,SPOUT(J,I,1),SPOUT(J,I,2),
     1 SPOUT(J,I,3),SPOUT(J,I,4),SPOUT(J,I,5)
          WRITE(12,933)SPOUT(J,I,6),SPOUT(J,I,7)

        ENDIF


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
