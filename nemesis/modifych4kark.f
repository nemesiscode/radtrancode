      subroutine modifych4kark(npro,patm,temp,ch4tropvmr,ch4stratvmr,
     1  RH,slope,xnew,xnewgrad)
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INTEGER npro,I,J,K,J3
      REAL PATM(MAXPRO),TEMP(MAXPRO)
      REAL XNEWGRAD(MAXPRO),dslope,slopet
      REAL ch4stratvmr,ch4tropvmr,RH,slope,pbar(maxpro),psvp(maxpro)
      REAL xref(maxpro),xnew(maxpro),slopex,dp,dch4,xnew1(maxpro)
      REAL tmp,pch4(maxpro)
      real SCH40,SCH41

C     Thermodynamic sata extracted from Handbook of Physics and Chemistry and
C     Kaye and Labyfor PRAXIS book. 
C     Taken from /home/irwinp/PRAXIS/properties/gravity/calc_wlapse.pro
C     SCH4=[10.6815,-1163.83] ; sublimation L'*R = 9671.42
      parameter (SCH40=10.6815,SCH41=-1163.83)      
C     NB, psvp is in bar.
      j3=-1

      print*,ch4tropvmr,ch4stratvmr,RH
      do 10 i = 1,npro
        pbar(i)=patm(i)*1.013
        if(pbar(i).gt.0.3)then
         j3=i
        endif
        tmp=SCH40+SCH41/temp(i)
        if(tmp.lt.-69.0)then
         psvp(i)=1e-30
        else
         psvp(i)=exp(tmp)
        endif

        pch4(i)=ch4tropvmr*pbar(i)      
        if(pch4(i)/psvp(i).gt.1.0)then
         pch4(i)=psvp(i)*RH
        endif
        if(pbar(i).lt.0.1.and.pch4(i)/pbar(i).gt.ch4stratvmr)then
         pch4(i)=pbar(i)*ch4stratvmr
        endif
        if(pbar(i).gt.0.5.and.pch4(i)/pbar(i).gt.ch4tropvmr)then
         pch4(i)=pbar(i)*ch4tropvmr
        endif
        xref(i)=pch4(i)/pbar(i)
        xnew(i)=xref(i)
        xnew1(i)=xref(i)
        print*,i,psvp(i),patm(i),pbar(i),xref(i)
10    continue

      j3=j3+1

      k=-1
      do 20 j=j3,1,-1  
       dp = pbar(j)-pbar(j+1)
       dch4 = xref(j)-xref(j+1)
       if(dch4.ge.0.0) then
        slopex = 0.01*dp/dch4
        print*,j,pbar(j),xref(j),xnew(j),dp,dch4,slopex,slope
        if(slopex.lt.slope.and.k.lt.0) then
         k=j
        endif
       endif
20    continue

25    print*,'kturn=',k
      print*,pbar(k),pbar(k+1),pbar(k)-pbar(k+1)
      print*,xref(k),xref(k+1),xref(k)-xref(k+1)
      print*,0.01*dp/dch4,slope

      if(k.ge.0) then
       do 30 j=k,1,-1
        dp = pbar(j)-pbar(j+1)
        xnew(j)= xnew(j+1)+0.01*dp/slope
        if(xnew(j).gt.ch4tropvmr)then
         xnew(j)=ch4tropvmr
        endif
        print*,j,dp,xref(j),xnew(j+1)+0.01*dp/slope,xnew(j)
30     continue
      endif


      do i=1,npro
       xnewgrad(i)=0.
      enddo

      dslope=0.1
      slopet = slope+dslope

      k=-1
      do 21 j=j3,1,-1
       dp = pbar(j)-pbar(j+1)
       dch4 = xref(j)-xref(j+1)
       if(dch4.gt.0.0) then
        slopex = 0.01*dp/dch4
        if(slopex.lt.slopet.and.k.lt.0) then
         k=j
        endif
       endif
21    continue

      if(k.ge.0) then
       do 31 j=k,1,-1
        dp = pbar(j)-pbar(j+1)
        xnew1(j)= xnew1(j+1)+0.01*dp/slopet
        if(xnew1(j).gt.ch4tropvmr)then
         xnew1(j)=ch4tropvmr
        endif
        xnewgrad(i)=(xnew1(i)-xnew(i))/dslope
31     continue
      endif

      return

      end
