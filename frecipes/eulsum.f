      SUBROUTINE EULSUM(SUM,TERM,JTERM,WKSP)
      DIMENSION WKSP(JTERM)
      IF(JTERM.EQ.1)THEN
        NTERM=1
        WKSP(1)=TERM
        SUM=0.5*TERM
      ELSE
        TMP=WKSP(1)
        WKSP(1)=TERM
        DO 11 J=1,NTERM
          DUM=WKSP(J+1)
          WKSP(J+1)=0.5*(WKSP(J)+TMP)
          TMP=DUM
11      CONTINUE
        IF(ABS(WKSP(NTERM+1)).LE.ABS(WKSP(NTERM)))THEN
          SUM=SUM+0.5*WKSP(NTERM+1)
          NTERM=NTERM+1
        ELSE
          SUM=SUM+WKSP(NTERM+1)
        ENDIF
      ENDIF
      RETURN
      END
