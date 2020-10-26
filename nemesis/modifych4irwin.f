      subroutine modifych4irwin(npro,patm,temp,ch4tropvmr,ch4stratvmr,
     1  RH,xnew,xnewgrad)
C     ************************************************************************************
C     Subroutine set methane profile as defined by Irwin et al. (2020)
C 
C     Input variables:
C	npro		integer	Number of levels in profile
C	patm(npro)	REAL	Pressure of atmosphere (atm)
C	temp(npro)	REAL	Pressure of atmosphere (atm)
C	ch4tropvmr	REAL	Required limiting tropospheric VMR
C	ch4stratvmr	REAL	Required limiting stratospheric VMR
C	RH		REAL	Required limiting relative humidity (fraction, i.e., in range 0 to 1)
C
C     Output variables
C	xnew(npro)	REAL	Output CH4 VMR profile
C	xnewgrad(npro)	REAL	Rate of change of CH4 VMRs with deep CH4 vmr, ch4tropvmr
C
C     Pat Irwin
C	20/6/20
C
C     ************************************************************************************

      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INTEGER npro,I,J,K
      REAL PATM(MAXPRO),TEMP(MAXPRO)
      REAL XNEWGRAD(MAXPRO)
      REAL ch4stratvmr,ch4tropvmr,RH,pbar(maxpro),psvp(maxpro)
      REAL xnew(maxpro)
      REAL tmp,pch4(maxpro)
      real SCH40,SCH41

C     Thermodynamic sata extracted from Handbook of Physics and Chemistry and
C     Kaye and Labyfor PRAXIS book. 
C     Taken from /home/irwinp/PRAXIS/properties/gravity/calc_wlapse.pro
C     SCH4=[10.6815,-1163.83] ; sublimation L'*R = 9671.42
      parameter (SCH40=10.6815,SCH41=-1163.83)      
C     NB, psvp is in bar.

      do 10 i = 1,npro
        xnewgrad(i)=0.
        pbar(i)=patm(i)*1.013
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
         xnewgrad(i)=1.0
        endif
        xnew(i)=pch4(i)/pbar(i)
10    continue

      return

      end
