      subroutine modifynh3irwin(npro,patm,temp,NH3tropvmr,NH3stratvmr,
     1  RH,xnew,xnewgrad)
C     ************************************************************************************
C     Subroutine set NH3 profile in the same was as the methane profile defined by Irwin et al. (2020)
C 
C     Input variables:
C	npro		integer	Number of levels in profile
C	patm(npro)	REAL	Pressure of atmosphere (atm)
C	temp(npro)	REAL	Pressure of atmosphere (atm)
C	NH3tropvmr	REAL	Required limiting tropospheric VMR
C	NH3stratvmr	REAL	Required limiting stratospheric VMR
C	RH		REAL	Required limiting relative humidity (fraction, i.e., in range 0 to 1)
C
C     Output variables
C	xnew(npro)	REAL	Output NH3 VMR profile
C	xnewgrad(npro)	REAL	Rate of change of NH3 VMRs with deep NH3 vmr, NH3tropvmr
C
C     Pat Irwin
C	20/6/20
C       11/7/23         Modified for NH3.
C
C     ************************************************************************************

      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INTEGER npro,I,J,K
      REAL PATM(MAXPRO),TEMP(MAXPRO)
      REAL XNEWGRAD(MAXPRO)
      REAL NH3stratvmr,NH3tropvmr,RH,pbar(maxpro),psvp(maxpro)
      REAL xnew(maxpro)
      REAL tmp,pNH3(maxpro)
      real SNH30,SNH31
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C     Thermodynamic sata extracted from Handbook of Physics and Chemistry and
C     Kaye and Labyfor PRAXIS book. 
C     Taken from /home/irwinp/PRAXIS/properties/gravity/calc_wlapse.pro
C     SNH3=[17.3471,-3930.55] ; sublimation L'*R = 32662.87
      parameter (SNH30=17.3471,SNH31=-3930.55)      
C     NB, psvp is in bar.

      do 10 i = 1,npro
        xnewgrad(i)=0.
        pbar(i)=patm(i)*1.013
        tmp=SNH30+SNH31/temp(i)
        if(tmp.lt.-69.0)then
         psvp(i)=1e-30
        else
         psvp(i)=exp(tmp)
        endif

        pNH3(i)=NH3tropvmr*pbar(i)      
        if(pNH3(i)/psvp(i).gt.1.0)then
         pNH3(i)=psvp(i)*RH
        endif
        if(pbar(i).lt.0.1.and.pNH3(i)/pbar(i).gt.NH3stratvmr)then
         pNH3(i)=pbar(i)*NH3stratvmr
        endif
        if(pbar(i).gt.0.5.and.pNH3(i)/pbar(i).gt.NH3tropvmr)then
         pNH3(i)=pbar(i)*NH3tropvmr
         xnewgrad(i)=1.0
        endif
        xnew(i)=pNH3(i)/pbar(i)
10    continue

      return

      end
