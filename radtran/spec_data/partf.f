      REAL FUNCTION PARTF(ID,ISO,TEMP,IPTF)
C     $Id: partf.f,v 1.3 2011-06-17 15:53:02 irwin Exp $
C-----------------------------------------------------------------------------
C_TITLE:  PARTF: computes total partition function
C
C_ARGS:   ID:INTEGER     local gas identifier
C         ISO:INTEGER    local isotope identifier
C         TEMP:REAL      temperature in Kelvin
C	  IPTF:INTEGER   Flag to indicate alternative processing
C
C_KEYS:   FUNCT,ATMOS,LBL
C
C_DESCR:  uses four term polynomial fits to compute ration of total partition
C         functions for use by line by line programs.
C         Representation is as described in Gimache et al 1990 and distributed
C         with HITRAN 91
C
C_DESCR1: The partition function is used to scale the line strengths to 
C	  different temperatures. It has two components, rotational and
C	  vibrational. Thus the line strength ratio is:
C
C	  S(T)      Z_v(T0)*Z_r(T0)
C	  ---    =  --------------  x other terms.
C         S(T0)      Z_v(T)*Z_r(T)
C
C         The rotational partition function is relatively easy to
C         calculate and the above may be simplified to:
C
C                       n
C	  S(T)      (T0)   Z_v(T0)
C	  ---    =  (--) x ------  x other terms.
C         S(T0)     ( T)   Z_v(T)
C
C         where n = 1.0 for linear molecules and 1.5 for non-linear
C         molecules
C
C         To fold in the vibrational part of the partition function
C         we use a 4th order polynomial scheme:
C
C	  S(T)      A + B*T0 + C*T0^2 + D*T0^3
C	  ---    =  --------------------------  x other terms.
C         S(T0)     A + B*T  + C*T^2  + D*T^3
C         
C         where A,B,C,D are coefficients which are read from the
C	  'gas_info.dat' file.
C
C         for linear molecules with no great vibrational contrib to the 
C	  partition function we set:
C
C		A=0
C		B=1
C		C=0
C		D=0
C	  
C         for linear molecules with no great vibrational contrib to the 
C         partition function we set:
C
C		A= -0.30755E+02 
C		B = 0.63608E+01 
C		C = 0.43890E-03
C		D = -0.23683E-06
C
C         which is a reasonable fit to T**1.5
C 
C_CALLS:  
C
C_BUGS:   not all gases have valid fits
C
C_HIST:  23oct92  SBC ORIGINAL VERSION with values copied from GENLN2 QTIPS
C        30JUN10  JML Updated CH4 coefficients for three temperature ranges
C-----------------------------------------------------------------------------
      INCLUDE '../includes/dbcom.f'
C-----------------------------------------------------------------------------
      INTEGER ID,ISO
      REAL TEMP
C
      REAL T0,T02,T03
      PARAMETER (T0=296.,T02=T0*T0,T03=T0*T02)
C
      REAL PARTT0,PTMP

c     ch4poly(temprange,iso,poly coeff)
      real ch4poly(3,3,4)
      data (ch4poly(1,1,i),i=1,4,1)/-3.14109e+02, 4.59494e+00,
     1 -8.39981e-03,1.19895e-05/
      data (ch4poly(1,2,i),i=1,4,1)/-6.28187e+02, 9.18943e+00,
     1 -1.67991e-02,2.39779e-05/
      data (ch4poly(1,3,i),i=1,4,1)/-2.54542e+03, 3.71886e+01,
     1 -6.79955e-02,9.70626e-05/
      data (ch4poly(2,1,i),i=1,4,1)/-2.60101e+05, 6.48085e+02,
     1 -5.45139e-01,1.63916e-04/
      data (ch4poly(2,2,i),i=1,4,1)/-5.20196e+05, 1.29615e+03,
     1 -1.09025e+00,3.27821e-04/
      data (ch4poly(2,3,i),i=1,4,1)/-2.10531e+06, 5.24595e+03,
     1 -4.41287e+00,1.32696e-03/
      data (ch4poly(3,1,i),i=1,4,1)/-1.19710e+07, 1.68635e+04,
     1 -8.07604e+00,1.33856e-03/
      data (ch4poly(3,2,i),i=1,4,1)/-2.39334e+07, 3.37157e+04,
     1 -1.61471e+01,2.67637e-03/
      data (ch4poly(3,3,i),i=1,4,1)/-9.68816e+07, 1.36480e+05,
     1 -6.53634e+01,1.08341e-02/


      PARTT0=DBQTA(ISO,ID)+DBQTB(ISO,ID)*T0
     1  +DBQTC(ISO,ID)*T02+DBQTD(ISO,ID)*T03
C
      PTMP=PARTT0/(DBQTA(ISO,ID) + TEMP*( DBQTB(ISO,ID) +
     1  TEMP*( DBQTC(ISO,ID) + TEMP*DBQTD(ISO,ID) ) ) )

      IF(IPTF.EQ.1.AND.ID.EQ.6)THEN 
       PARTT0=ch4poly(1,ISO,1)+ch4poly(1,ISO,2)*T0
     1+ch4poly(1,ISO,3)*T02+ch4poly(1,ISO,4)*T03

	if (TEMP.LE.1000.0) then
      PTMP=PARTT0/(ch4poly(1,ISO,1) + TEMP*( ch4poly(1,ISO,2) +
     1TEMP*( ch4poly(1,ISO,3) + TEMP*ch4poly(1,ISO,4) ) ) )	

	else if (TEMP.GT.1000.0.AND.TEMP.LE.2000.0) then
      PTMP=PARTT0/(ch4poly(2,ISO,1) + TEMP*( ch4poly(2,ISO,2) +
     1TEMP*( ch4poly(2,ISO,3) + TEMP*ch4poly(2,ISO,4) ) ) )

	else if (TEMP.GT.2000) then
      PTMP=PARTT0/(ch4poly(3,ISO,1) + TEMP*( ch4poly(3,ISO,2) +
     1TEMP*( ch4poly(3,ISO,3) + TEMP*ch4poly(3,ISO,4) ) ) )

	end if

      ENDIF

      IF(PTMP.LT.0.0)THEN
       PRINT*,'PROBLEM IN PARTF - RESULT IS NEGATIVE'
       PRINT*,'ID,ISO,IPTF = ',ID,ISO,IPTF
       PRINT*,'TEMP = ',TEMP
       PRINT*,'DBQTA,DBQTB,DBQTC,DBQTD',DBQTA(ISO,ID),DBQTB(ISO,ID),
     &DBQTC(ISO,ID),DBQTD(ISO,ID)
       PRINT*,'PARTT0 = ',PARTT0
       PRINT*,'PARTF = ',PTMP
       STOP
      ENDIF

      PARTF=PTMP

      RETURN
      END
