      subroutine map2pro(nparam,doutputdq,doutmoddq)
      IMPLICIT NONE
C     *****************************************************************
C     Subroutine which converts from rate of change of output with respect
C     to layer properties to rate of change of output with .prf properties
C
C     Input variables
C	nparam		integer		Number of variables
C	doutputdq(maxpath,maxlay,maxgas+2+maxcon) real Rate of change
C			 		of output with respect to layer
C					properties:
C
C	1 to ngas			layer gas amounts (mol cm-2)
C	ngas+1				layer temperatures
C	ngas+2 to ngas+2+ncont-1	layer cloud densities
C	ngas+ncont+2			layer para-H2 fraction.
C
C     Output variables
C	doutmoddq(maxpat,maxgas+2+maxcon,maxpro) real Rate of change
C					of output with respect to .prf
C					profile properties:
C
C	1 to ngas			profile gas vmrs
C	ngas+1				profile temperatures
C	ngas+2 to ngas+2+ncont-1	profile cloud specific densities
C	ngas+ncont+2			profile para-H2 fraction.
C
C     The routine multipies the gradients calculated per layer by 
C     the matrices DTE, DFP, DCO and DAM which convert the gradients with
C     respect to temperature, para-H2 fraction, dust opacity and layer
C     amount to the profile temperatures, para-H2 fraction, specific
C     concentration and vmr respectively.
C
C     Pat Irwin		30/7/01
C     Pat Irwin		5/12/03	Removed references to DEMTEMP. Use DTE
C				instead.
C     Pat Irwin		29/2/12	Updated for Radtrans2.0
C	
C     *****************************************************************

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/laycom.f'
      INCLUDE '../includes/laygrad.f'
      REAL doutputdq(maxpat,maxlay,maxgas+2+maxcon)
      REAL doutmoddq(maxpat,maxgas+2+maxcon,maxpro)
      INTEGER nparam,ipath,k,jj,j,i,nlay1,ipath1

      do ipath=1,npath
        do k=1,nparam
          do jj=1,npro
            doutmoddq(ipath,k,jj)=0.0
C           Find correct number of atmospheric layers for IMOD=8 
C		(i.e. SCR radiometer). 
            nlay1=nlayin(ipath)
            ipath1=ipath
            if(imod(ipath).eq.8)then
              nlay1=nlayin(1)
              ipath1=1
            endif
            do j=1,nlay1
             if(k.le.ngas)then  
              doutmoddq(ipath,k,jj)=doutmoddq(ipath,k,jj)+
     1          doutputdq(ipath,j,k)*DAM(LAYINC(J,IPATH1),JJ)
             endif

             if(k.eq.ngas+1)then
              doutmoddq(ipath,k,jj)=doutmoddq(ipath,k,jj)+
     1          doutputdq(ipath,j,k)*DTE(LAYINC(J,IPATH1),JJ)
             endif

             if(k.gt.ngas+1)then
              if((k-ngas-1).le.ncont) then
                doutmoddq(ipath,k,jj)=doutmoddq(ipath,k,jj)+
     1            doutputdq(ipath,j,k)*DCO(LAYINC(J,IPATH1),JJ)
              else
                doutmoddq(ipath,k,jj)=doutmoddq(ipath,k,jj)+
     1            doutputdq(ipath,j,k)*DFP(LAYINC(J,IPATH1),JJ)
              endif
             endif
            enddo
          enddo
        enddo
      enddo

      return

      end
