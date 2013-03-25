      SUBROUTINE INTERP_PT(NPRO,H,P,T,HEIGHT,PNOW,TNOW,F,IFL)
C     **************************************************************
C     Interpolate pressure and temperature profiles to find current
C     pressure, temperature and position in profile arrays in terms
C     of position in array and fraction.
C
C     Input variables
C  	NPRO	INTEGER	Number of vertical levels in profile
C	H(NPRO)	REAL	Profile heights
C	P(NPRO)	REAL	Profile pressures
C	T(NPRO)	REAL	Profile temperatures
C	HEIGHT	REAL 	Current height
C
C     Output variables
C	PNOW	REAL	Current pressure
C	TNOW	REAL	Current temperature
C	F	REAL	Fractional position in array relative to nearest
C			value below
C	IFL	REAL	Nearest value below. 
C			i.e. HEIGHT= (1-F)*H(IFL)+F*H(IFL+1)
C
C     Pat Irwin	22/3/13
C
C     **************************************************************
      IMPLICIT NONE
      INTEGER NPRO,IFL,K
      REAL HEIGHT,H(NPRO),PNOW,P(NPRO),TNOW,T(NPRO),F
      F=-1.0
      IFL = 0
      DO K=1,NPRO-1
        IF(HEIGHT.GE.H(K).AND.HEIGHT.LT.H(K+1))THEN
         F = (HEIGHT - H(K))/(H(K+1)-H(K))
         IFL=K
        ENDIF
      ENDDO
      IF(IFL.EQ.0)THEN
        IF(HEIGHT.LT.H(1))THEN
         IFL=1
         F=0.0
        ENDIF
        IF(HEIGHT.GE.H(NPRO))THEN
         IFL=NPRO-1
         F=1.0
        ENDIF
      ENDIF


      PNOW = (1-F)*P(IFL) + F*P(IFL+1)
      TNOW = (1-F)*T(IFL) + F*T(IFL+1)


      RETURN

      END

