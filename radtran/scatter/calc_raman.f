      subroutine calc_raman(v,T,inormal,xrot,xvib)
C     ******************************************************************
C     Subroutine to determine effective Raman scattering opacities in the 
C     approximation that the temperature is cold and thus H2 exists only in
C     rotational states J=0 and J=1. 
C     
C     Taken from Sromovsky (2005) Icarus 173, 254-283, Eq. 3
C
C     Input variables
C	v	real	Wavenumber (cm-1)
C	T	real	Temperature (K)
C	inormal	integer	0 = equilibrium ortho/para, 1= 3:1 ortho/para
C
C     Output 
C	xrot	real	Rotational Raman x-section (cm2)
C	xvib	real	Vibrational Raman x-section (cm2)
C
C     Pat Irwin	2/12/16	Original version
C 
C     ******************************************************************

      real v,xrot,xvib,Qout(20),T,fpara
      double precision T1

      call ramanh2(v,Qout)
      T1=dble(T)
      if(inormal.eq.0)then
       fpara=calcpara(T1)
      else
       fpara=0.25
      endif
      xrot=fpara*Qout(2)+(1.0-fpara)*Qout(6)
      xvib=0.5*(Qout(3)+Qout(7)+Qout(4)+Qout(8))

      return

      end
   
