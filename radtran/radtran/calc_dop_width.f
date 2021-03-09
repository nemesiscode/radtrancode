      REAL FUNCTION CALC_DOP_WIDTH(ABSCO,AD,DOP_MOD,RODGERS)
C     $Id:
C     ***********************************************************
C     Calculates the Doppler width of a line assuming different
C     assumptions.
C
C     Pat Irwin	?/?/??	Original
C
C     ***********************************************************
      REAL ABSCO,AD,PI
      PARAMETER(PI=3.1415927)
      DOUBLE PRECISION RODGERS(4,0:7),SUM,Z
      INTEGER DOP_MOD
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


      IF(DOP_MOD.EQ.0)THEN
C     Weak limit
       CALC_DOP_WIDTH = ABSCO

      ELSEIF(DOP_MOD.EQ.1)THEN
C     Strong Doppler equivalent width
       CALC_DOP_WIDTH = 2*AD*(LOG(0.5641895*ABSCO/AD))**0.5

      ELSEIF(DOP_MOD.EQ.2)THEN
C     Combined Doppler width after Rodgers and Williams 1974
       Z=DBLE(ABSCO/(AD*SQRT(PI)))
       IF(Z.LE.0.)THEN
        CALC_DOP_WIDTH=0.
        if(idiag.gt.0)print*,'Error. Calc_Dop_Width z.le.0'
       ELSE
        SUM=0.
        IF(Z.LE.5.)THEN
         DO 10 N=0,7
          SUM=SUM+RODGERS(3,N)*Z**N
10       CONTINUE
         SUM=SUM*Z*SQRT(PI)
        ELSE
         DO 20 N=0,7
          SUM=SUM+RODGERS(4,N)*(DLOG(Z)**(-N))
20       CONTINUE
         SUM=SUM*DSQRT(DLOG(Z))
        ENDIF
        CALC_DOP_WIDTH=AD*SNGL(SUM)
       END IF
      END IF

      RETURN

      END


