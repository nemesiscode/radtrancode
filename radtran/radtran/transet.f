      subroutine transet(bandtyp,tpout,jgas,nav,iav,w,qrot,tpart,
     1 P,T,q,x,y,sig,ndata)
C     ********************************************************************
C     Subroutine to compute a transmission function (T as a function of 
C     absorber amount) for band data averaged over a specified interval.
C
C     Extracted from esumseta.f to make code more modular
C
C     input variables
C	bandtyp(mgas)		integer Band model type
C	tpout(mbin,mgas,7)	real	Band data
C	jgas			integer	Which gas id (in tpout) to use
C	nav			integer	Number of bins to average
C	iav(mbin)		integer	Which bins
C	w(mbin)			real	Weight of bin in final average
C	qrot			real	Rotation Z T-exponent
C	P			real	Pressure (atm)
C	T			real	Temperature (K)
C	q			real	mole fraction (0=fb,1=sb)
C	ndata			integer	Number of points in required
C					transmission function
C	imod			integer Band model flag:
C					 0 = Kim Strong
C					 1 = Pat Irwin (corrected V)
C                                        2 = ???
C                                        3 = Karkoschka formulation
C
C     Output variables
C	x(ndata)		real	Absorber amounts
C	y(ndata)		real	Transmissions
C	sig(ndata)		real	Required fitting error
C
C     Pat Irwin		23/1/04
C
C     ********************************************************************

      implicit none
      integer mbin,mgas,jgas,i1,k,mdata
      parameter(mbin=8000,mgas=20,mdata=20)
      real qrot,P,T,q,x(mdata),y(mdata)
      real sig(mdata),dpexp,w(mbin),tpart(4)
      real tpout(mbin,mgas,7),C1,C2
      integer bandtyp(mgas),j
      real U,T1,U1,U2,DU,E_TRAN,yv
      real pi
      parameter (pi=3.1415927)
      integer NDATA,I,ng,imod,icalc,j05,iav(mbin),nav
      double precision SL,BL


     

C      print*,'TRANSET called'
C      print*,jgas,nav
C      print*,'bandtype : ',bandtyp(jgas)
C      do i=1,nav
C       print*,i,iav(i),w(i)
C      enddo
C      print*,qrot,P,T,q
C      print*,'tpart : ',(tpart(i),i=1,4)
C      print*,ndata

C     *********** first find where Transmission approximately 0.5
      U=1.

122   call calc_tau(bandtyp,tpout,jgas,nav,iav,w,U,P,T,q,qrot,
     1  tpart,T1)

      PRINT*,'U,T1',U,T1
      IF(T1.LT.0.5)THEN
        U=U/2.
        GOTO 122
      END IF

123   call calc_tau(bandtyp,tpout,jgas,nav,iav,w,U,P,T,q,qrot,
     1  tpart,T1)

      print*,'U,T1',U,T1
      IF(T1.GT.0.5)THEN
        IF(U.LE.1.0E38)THEN
         U=U*2.
         GOTO 123
        ENDIF
      END IF

C     ************ find upper limit of curve
      U2=U
      print*,'U for T=0.5 = ',U
222   call calc_tau(bandtyp,tpout,jgas,nav,iav,w,U2,P,T,q,qrot,
     1  tpart,T1)

      IF(T1.GT.0.001.AND.U2.LT.1E36)THEN
        IF(U2.LE.1.0E38)THEN
         U2=U2*2.
         GOTO 222
        ENDIF
      END IF

      print*,'U2,T1',U2,T1

C     *********** find lower limit of curve
      U1=U
322   call calc_tau(bandtyp,tpout,jgas,nav,iav,w,U1,P,T,q,qrot,
     1  tpart,T1)
      IF(T1.LT.0.999.AND.U1.GT.1E-36)THEN
       	 U1=U1/2.
       	 GOTO 322
      END IF

      PRINT*,'U1,T1',U1,T1

C     *********** Evaluate curve
      U1=LOG(U1)
      U2=LOG(U2)
      DU=(U2-U1)/(NDATA-1)
      E_TRAN=0.001

      DO 101 I=1,NDATA
       	 U=U1+DU*(I-1)
       	 U=DPEXP(U)
         call calc_tau(bandtyp,tpout,jgas,nav,iav,w,U,P,T,q,qrot,
     1    tpart,T1)
	 Y(I) = T1
       	 X(I)=U
       	 SIG(I)=E_TRAN
         print*,I,X(I),Y(I),SIG(I)
101   CONTINUE

      RETURN

      END

 
      subroutine calc_tau(bandtyp,tpout,jgas,nav,iav,w,U,P,T,q,
     1 qrot,tpart,T1)
C     ****************************************************************
C     Subroutine to calculate band-averaged transmission over a number 
C     of bins of band data
C
C     Input variables
C	bandtyp(mgas)		integer Band model type:
C					 -1 = T_GV - EKS
C					 0 = T_GV - extended voigt
C					 1 = T_GV - Kam C1,C2
C					 2 = T_GV, 2-EL, T**QROT
C					 3 = T_GV, 2-EL, tpart(4)
C					 4 = Karkoschka09
C					 5 = Goody-Lorentz
C	tpout(mbin,mgas,7)	real	Band data
C	jgas			integer	Which gas id (in tpout) to use
C	nav			integer	Number of bins to average
C	iav(mbin)		integer	Which bins
C	w(mbin)			real	Weight of bin in final average
C	U			real	absorber amount
C	P			real	Pressure (atm)
C	T			real	Temperature (K)
C	q			real	mole fraction (0=fb,1=sb)
C	qrot			real	Rotation Z T-exponent
C	tpart(4)		real	Partition function coeffs (for
C					imod=1)
C	imod			integer Band model flag:
C					 0 = Pat Irwin (corrected V)
C     					 1 = As 0 but using tpart
C       icalc			integer Band model calculation ID
C					 0 = single EL
C					 1  = a, E1, E2
C     Output variables
C	T1			real	Transmission
C
C     Pat Irwin		23/1/04
C     ****************************************************************
      implicit none
      integer i,jgas,mgas,mbin,imod,i1
      parameter (mgas=20,mbin=8000)
      real knu0,delad,y0,EL,E1,E2,A,SFB,qrot,tpart(4),P,T,q
      real T1,xcorr,U,tautmp,dpexp,w(mbin),C1,C2,SUMT
      integer NT,iav(mbin),nav,bandtyp(mgas),j
      real TGV3,TAU_EKS,TKARK,YV,TAU_MG_LOR
      real tpout(mbin,mgas,7),ALCORR,TX1,TX2,WX1,WX2
      real p0,t0,kap100,kap198,kap296,dline
      parameter (p0=1.,t0=296.)
      logical MG

      T1 = 0.0
      SUMT = 0.0
      do 101 i1=1,nav
        I = IAV(I1)
        IF(BANDTYP(JGAS).LT.4)THEN
          KNU0 = TPOUT(I,JGAS,1)
          DELAD = TPOUT(I,JGAS,2)
          Y0 = TPOUT(I,JGAS,3)
          E1 = TPOUT(I,JGAS,4)
          SFB = TPOUT(I,JGAS,5)
          E2 = TPOUT(I,JGAS,6)
          A = TPOUT(I,JGAS,7)
C          print*,KNU0,DELAD,Y0,E1,SFB
          IF(KNU0.EQ.0)THEN
           tautmp = 0.0
          ELSE
           XCORR=0.0
           IF(BANDTYP(JGAS).LT.2)THEN
            IMOD=1
            A=0.5
            E2=E1
            IF(BANDTYP(JGAS).EQ.1)THEN
             C1 = ABS(TPOUT(I,J,6)*1E-4)
             C2 = TPOUT(I,J,7)
             XCORR = C1*U*(P/P0)*(T/T0)**C2
            ENDIF
           ELSE
            IF(BANDTYP(JGAS).EQ.2)THEN
             IMOD=1
            ELSE
             IMOD=2
            ENDIF
           ENDIF
           IF(BANDTYP(JGAS).LT.0)THEN
             C1 = ABS(TPOUT(I,J,6)*1E-4)
             C2 = TPOUT(I,J,7)
             XCORR = C1*U*(P/P0)*(T/T0)**C2
             tautmp = XCORR + TAU_EKS(KNU0,DELAD,Y0,E1,SFB,QROT,P,T,
     1          U,Q)
           ELSE
             tautmp = XCORR + TGV3(IMOD,QROT,TPART,KNU0,DELAD,Y0,A,
     1          E1,E2,SFB,P,T,U,Q)
C             print*,XCORR,tautmp
           ENDIF
          ENDIF

        ELSE

          IF(BANDTYP(JGAS).EQ.4)THEN

           KAP100 = TPOUT(I,JGAS,1)
           KAP198 = TPOUT(I,JGAS,2)
           KAP296 = TPOUT(I,JGAS,3)
           DLINE  = TPOUT(I,JGAS,4)
C          print*,kap100,kap198,kap296,dline
C          print*,P,T,U,Q
          
           tautmp = TKARK(KAP100,KAP198,KAP296,
     1      DLINE,P,T,U,Q)
C          print*,'tautmp = ',tautmp
          ELSE
           
           KNU0 = TPOUT(I,JGAS,1)
           DELAD = TPOUT(I,JGAS,2)
           Y0 = TPOUT(I,JGAS,3)
           EL = TPOUT(I,JGAS,4)
           SFB = TPOUT(I,JGAS,5)

           MG=.FALSE.
           ALCORR=1.
 
           YV = ALCORR*Y0/DELAD

           tautmp = TAU_MG_LOR(MG,KNU0,YV,EL,SFB,QROT,P,T,U,Q)
C          print*,'knu0,yv,el,sfb,qrot,tautmp = ',knu0,yv,el,sfb,
C     1     qrot,tautmp

          ENDIF
        ENDIF

        T1 = T1+w(i1)*DPEXP(-tautmp)
        SUMT=SUMT+w(i1)
C        print*,'sum',i1,w(i1),tautmp,T1,SUMT

101     continue

        T1=T1/SUMT

        if(T1.LT.0.0.OR.T1.GT.1.0)THEN
         print*,'calc_tau: transmission out of range'
         print*,'Transmission = ',T1
         print*,'bandtyp = ',bandtyp
         print*,'jgas,nav = ',jgas,nav
         print*,'XCORR etc. ',XCORR,IMOD,QROT,TPART,KNU0,DELAD,Y0,A,
     1          E1,E2,SFB,P,T,U,Q
         print*,tautmp,TGV3(IMOD,QROT,TPART,KNU0,DELAD,Y0,A,
     1          E1,E2,SFB,P,T,U,Q)
         do i=1,nav
          print*,iav(i),w(i)
         enddo
         print*,'tautmp,dpexp(-tautmp),knu0',tautmp,dpexp(-tautmp),knu0
         print*,'U,P,T,q,qrot,tpart',U,P,T,q,qrot,tpart
         print*,'Normalisation = ',SUMT
         if(T1.gt.1.0) T1=1.0
         if(T1.lt.0.0) T1=0.0
         stop
        ENDIF

        return

        end
