      subroutine ramanh2(v,Qout)
C     ***************************************************************
C     Subroutine to return the Raman scattering cross sections for molecular
C     hydrogen, computed using the fits of 
C     Ford and Browne, Atomic Data, 5, 305-313, 1973
C
C     Transitions are:                  Srom_ID
C	1 : J=0, Pure Rayleigh          R(0) 
C	2 : J=0, Rotational S-branch    S(0)
C	3 : J=0, Vibrational Q-branch   Q(0)
C	4 : J=0, Vibrational S-Branch   S_1(0)
C
C	5 : J=1, Pure Rayleigh          R(1)
C	6 : J=1, Rotational S-branch    S(1)
C	7 : J=1, Vibrational Q-branch   Q(1)
C	8 : J=1, Vibrational S-Branch   S_1(1)
C
C	9 : J=2, Pure Rayleigh
C	10: J=2, Rotational Q-branch
C	11: J=2, Rotational S-branch
C	12: J=2, Vibrational Q-branch
C	13: J=2, Vibrational O-branch
C	14: J=2, Vibrational S-Branch
C
C	15: J=3, Pure Rayleigh
C	16: J=3, Rotational Q-branch
C	17: J=3, Rotational S-branch
C	18: J=3, Vibrational Q-branch
C	19: J=3, Vibrational O-branch
C	20: J=3, Vibrational S-Branch
C
C     Input variables
C	v	real	Input wavenumber (cm-2)
C
C     Output variables
C	Qout(20) real	Computed cross-sections (cm2)
C
C     Pat Irwin
C     1/12/16

C     ***************************************************************
      implicit none
      integer ntype,itype,nc,I
      parameter(ntype=20,nc=9)
      double precision A(ntype,nc),Etrans(ntype),wl
      real v,Qout(ntype)
      double precision Q(ntype),xfact

      data (A(1,I),I=1,nc)/8.311d-45,0.0,1.244d-38,0.0,1.512d-32,
     &  0.0,1.724d-26,0.0,1.918d-20/
      data (A(2,I),I=1,nc)/2.563d-46,-2.119d-45,6.257d-40,-7.232d-39,
     &  1.079d-33,-1.619d-32,1.612d-27,-2.995d-26,2.228d-21/
      data (A(3,I),I=1,nc)/1.558d-46,-1.563d-44,4.588d-40,-6.486d-38,
     &  9.161d-34,-1.698d-31,1.551d-27,-3.585d-25,2.403d-21/
      data (A(4,I),I=1,nc)/3.335d-47,-5.359d-45,1.135d-40,-2.372d-38,
     &  2.592d-34,-6.676d-32,4.967d-28,-1.519d-25,8.610d-22/

      data (A(5,I),I=1,nc)/8.442d-45,0.0,1.275d-38,0.0,1.564d-32,0.0,
     &  1.801d-26,0.0,2.022d-20/
      data (A(6,I),I=1,nc)/1.550d-46,-2.154d-45,3.782d-40,-7.337d-39,
     &  6.523d-34,-1.640d-32,9.745d-28,-3.032d-26,1.348d-21/
      data (A(7,I),I=1,nc)/1.659d-46,-1.722d-44,4.957d-40,-7.217d-38,
     &  1.004d-33,-1.910d-31,1.724d-27,-4.076d-25,2.708d-21/
      data (A(8,I),I=1,nc)/2.033d-47,-3.427d-45,6.915d-41,-1.518d-38,
     &  1.582d-34,-4.280d-32,3.039d-28,-9.753d-26,5.282d-22/

      data (A(9,I),I=1,nc)/8.463d-45,0.0,1.279d-38,0.0,1.569d-32,0.0,
     & 1.805d-26,0.0,2.025d-20/
      data (A(10,I),I=1,nc)/5.105d-47,4.206d-46,1.242d-40,1.431d-39,
     & 2.136d-34,3.191d-33,3.180d-28,5.884d-27,4.380d-22/
      data (A(11,I),I=1,nc)/1.375d-46,-2.690d-45,3.399d-40,-9.274d-39,
     & 5.931d-34,-2.097d-32,8.961d-28,-3.917d-26,1.252d-21/
      data (A(12,I),I=1,nc)/1.653d-46,-1.700d-44,4.927d-40,-7.112d-38,
     & 9.955d-34,-1.880d-31,1.706d-27,-4.005d-25,2.674d-21/
      data (A(13,I),I=1,nc)/3.271d-48,-5.267d-46,1.302d-41,-2.543d-39,
     & 3.258d-35,-7.578d-33,6.619d-29,-1.796d-26,1.194d-22/
      data (A(14,I),I=1,nc)/5.604d-48,-1.317d-45,2.496d-41,-6.843d-39,
     & 6.761d-35,-2.157d-32,1.461d-28,-5.349d-26,2.775d-22/

      data (A(15,I),I=1,nc)/8.524d-45,0.0,1.295d-38,0.0,1.594d-32,0.0,
     & 1.840d-26,0.0,2.071d-20/
      data (A(16,I),I=1,nc)/6.622d-47,9.175d-46,1.611d-40,3.114d-39,
     & 2.770d-34,6.937d-33,4.123d-28,1.277d-26,5.680d-22/
      data (A(17,I),I=1,nc)/1.297d-46,-3.235d-45,3.211d-40,-1.116d-38,
     & 5.610d-34,-2.525d-32,8.488d-28,-4.723d-26,1.188d-21/
      data (A(18,I),I=1,nc)/1.667d-46,-1.714d-44,4.986d-40,-7.192d-38,
     & 1.011d-33,-1.906d-31,1.737d-27,-4.072d-25,2.730d-21/
      data (A(19,I),I=1,nc)/4.261d-48,-6.444d-46,1.699d-41,-3.111d-39,
     & 4.253d-35,-9.269d-33,8.634d-29,-2.195d-26,1.556d-22/
      data (A(20,I),I=1,nc)/3.947d-48,-1.057d-45,1.915d-41,-5.781d-39,
     & 5.475d-35,-1.887d-32,1.228d-28,-4.805d-26,2.401d-22/

      data (Etrans(I),I=1,8) /0.0,354.39,4162.06,4498.75,0.0,587.07,
     &  4156.15,4713.83/
      data (Etrans(I),I=9,14) /0.0,-354.39,814.48,4144.36,3807.67,
     & 4917.94/
      data (Etrans(I),I=15,20)/0.0,-587.07,1034.74,4126.76,3569.08,
     & 5109.35/

C     Find current wavelength in Angstroms
      wl=1e8/v

C     Now find different cross-sections (cm2)
      do 10 itype=1,ntype
       Q(itype)=0.
       do 20 I=1,nc
        Q(itype)=Q(itype)+A(itype,I)/wl**(I-1)
20     continue
       xfact=(1d8/wl-etrans(itype))**4
       if(xfact.gt.0.0)then
        Q(itype)=Q(itype)*xfact
       else
        Q(itype)=0.
       endif
       Qout(itype)=SNGL(Q(itype))
10    continue

      return

      end
