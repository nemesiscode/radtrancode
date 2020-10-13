      subroutine modifych4sromovsky(npro,patm,temp,ch4tropvmr,
     1  ch4stratvmr,PD,PT,RHC,RHM,VX,xnew)
C     ************************************************************************************
C     Subroutine set methane profile to descended methane profile of Sromovsky et al. (2019)
C 
C     Input variables:
C	npro		integer	Number of levels in profile
C	patm(npro)	REAL	Pressure of atmosphere (atm)
C	temp(npro)	REAL	Pressure of atmosphere (atm)
C	ch4tropvmr	REAL	Required limiting tropospheric VMR
C	ch4stratvmr	REAL	Required limiting stratospheric VMR
C	PD		REAL	Deep pressure (bar)
C	PT		REAL	Tropopause pressure (bar)
C	RHC		REAL	Cloud-top relative humidity fraction (0.0-1.0)
C	RHM		REAL	Tropopause minumum relative humidity (0.0-1.0)
C	VX		REAL	Descending parameter
C
C     Output variables
C	xnew(npro)	REAL	Output CH4 VMR profile
C
C     Pat Irwin
C	9/10/20
C
C     ************************************************************************************

      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INTEGER npro,I,J,K,ideep,icond,ipause,npro1
      REAL PATM(MAXPRO),TEMP(MAXPRO),LPBAR(MAXPRO),LP1(MAXPRO)
      REAL T1(MAXPRO),temp1(MAXPRO),ptmp(MAXPRO),alpha(MAXPRO)
      REAL Y(MAXPRO),alpha1(MAXPRO)
      REAL pfinebar(MAXPRO),tfine(MAXPRO),PC,PD,PT,RHC,RHM,VX
      REAL ch4stratvmr,ch4tropvmr,RH,pbar(maxpro),psvp(maxpro)
      REAL xnew(maxpro),a1,a2
      REAL tmp,pch4(maxpro),PM
      real SCH40,SCH41,SCH42

C     Thermodynamic data extracted from Handbook of Physics and Chemistry and
C     Kaye and Laby for PRAXIS book. 
C     Taken from /home/irwinp/PRAXIS/properties/gravity/calc_wlapse.pro
C     SCH4=[10.6815,-1163.83] ; sublimation L'*R = 9671.42
      parameter (SCH40=10.6815,SCH41=-1163.83)      
C     NB, psvp is in bar.

C     Data in SVP.DAT file
C     SCH4=[10.50377,-1116.966,-4.633456e-03] - need to check this!!
C      parameter (SCH40=10.50377,SCH41=-1116.966,SCH42=-4.633456e-03)      


C      print*,ch4tropvmr,ch4stratvmr,PD,PT,RHC,RHM,VX
      do i=1,npro
       j = npro-i+1
       pbar(i)=patm(i)*1.013
       lpbar(j)=alog(pbar(i))
       temp1(j)=temp(i)
C       print*,i,j,pbar(i),lpbar(j),temp1(j)
      enddo


C     First interpolate the TP profile onto a finer grid    
      npro1=50
      do 10 i=1,npro1
       LP1(I)=LPBAR(1)+(LPBAR(NPRO)-LPBAR(1))*FLOAT(I-1)/FLOAT(NPRO1-1)
       CALL VERINT(LPBAR,TEMP1,NPRO,T1(I),LP1(I))
       j = npro1-i+1
       pfinebar(j) = exp(lp1(i))
       tfine(j)=t1(i)
       ptmp(j)=pfinebar(j)
C       print*,i,j,pfinebar(j),lp1(i),tfine(j)
10    continue

C     First find level where CH4 saturates and limit CH4 profile to 100% RH
      icond=-1
      ipause=-1
      ideep=-1
      do 15 i = 1,npro1
        tmp=SCH40+SCH41/tfine(i)
C        tmp=SCH40+SCH41/tfine(i)+SCH42*tfine(i)
        if(tmp.lt.-69.0)then
         psvp(i)=1e-30
        else
         psvp(i)=exp(tmp)
        endif

        pch4(i)=ch4tropvmr*pfinebar(i)      
        if(pch4(i)/psvp(i).ge.1.0)then
         pch4(i)=psvp(i)
         if(icond.lt.0)icond = i
        endif
        if(pfinebar(i).le.PT.and.ipause.lt.0)then
         ipause=i
        endif
        if(pfinebar(i).le.PD.and.ideep.lt.0)then
         ideep=i
        endif
C        print*,'B',i,pfinebar(i),tfine(i),pch4(i),psvp(i),
C     &    pch4(i)/pfinebar(i)
15    continue

C      print*,'Meth. con. lev., press(bar), temp(K) = ',
C     1  icond,pfinebar(icond),tfine(icond)

      PC=pfinebar(icond)

C      print*,'ideep, icond, ipause = ',ideep, icond, ipause
C      print*,pfinebar(ideep),pfinebar(icond),pfinebar(ipause)

      PM=pfinebar(ipause)

C     Now need to apply Sromovsky relative humidity profile to region between cloud top and tropopause.
      a1=log(PC/PM)
      do 20 i=icond,ipause
        a2=log(PC/pfinebar(i))
        RH=RHM+(RHC-RHM)*(1-a2/a1)
C        print*,'D',i,a1,a2,RH,pch4(i),psvp(i)*RH
        pch4(i)=psvp(i)*RH
20    continue
      do 30 i=ipause+1,npro1
        pch4(i)=psvp(i)*RHM
        if(pch4(i)/pfinebar(i).gt.ch4stratvmr)then
         pch4(i)=pfinebar(i)*ch4stratvmr
        endif
30    continue

C      do i=1,npro1
C        print*,'C',i,pfinebar(i),tfine(i),pch4(i),psvp(i),
C     &    pch4(i)/pfinebar(i),pch4(i)/psvp(i)
C      enddo

C      print*,'A'
      do 35 i=1,npro1
       alpha(i)=pch4(i)/pfinebar(i)
C       print*,i,pfinebar(i),pch4(i),alpha(i)
35    continue

C     Now do descended profile correction
      do i=1,ideep-1
       ptmp(i)=pfinebar(1)
      enddo

      do 40 i=ideep,ipause
       a1 = (alpha(i)/ch4tropvmr)
       a2 = a1**VX
       ptmp(i)= pfinebar(i)*(1 + a2*(PD/PC-1))
       if(ptmp(i).gt.pfinebar(1))ptmp(i)=pfinebar(i)
C       print*,i,pfinebar(i),alpha(i),a1,a2,ptmp(i)
40    continue


C      do i=1,npro1
C       print*,i,pfinebar(i),ptmp(i),alpha(i)
C      enddo

      do 45 i=1,npro1
       j=npro1-i+1
       lp1(j)=log(ptmp(i))
       y(j)=alpha(i)
45    continue

      do 50 i=1,npro1
       call verint(lp1,y,npro1,alpha1(i),log(pfinebar(i)))
50    continue

      do 55 i=1,npro
       call verint(lp1,y,npro1,xnew(i),log(pbar(i)))   ! might need to reverse this!
55    continue

C      do 55 i=1,npro
C       xnew(i)=alpha(i)
C55    continue

      return
      
      end
