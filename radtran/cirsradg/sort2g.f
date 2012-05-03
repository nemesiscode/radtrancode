      SUBROUTINE SORT2G(N,RA,IB)
C     ******************************************************************
C     Modified numerical recipes routine to sort a vector RA of length N
C     into ascending order. Integer vector IB is initially set to
C     1,2,3,... and on output keeps a record of how RA has been sorted.
C
C     Pat Irwin	31/7/01	Original version
C     Pat Irwin	29/2/12	Updated for Radtrans2.0
C 
C     ******************************************************************
      include '../includes/arrdef.f'
      REAL RA(MAXRANK)
      INTEGER IB(MAXRANK)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          IRB=IB(L)
        ELSE
          RRA=RA(IR)
          IRB=IB(IR)
          RA(IR)=RA(1)
          IB(IR)=IB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            IB(1)=IRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            IB(I)=IB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        IB(I)=IRB
      GO TO 10
      END
