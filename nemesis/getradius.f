      REAL FUNCTION GETRADIUS(ICLOUD,NVAR,VARIDENT,VARPARAM,XN,NPRO)
C     **************************************************************
C     Function to sift through XN and VARIDENT FILES to find the mean 
C     radius of particles defined though IMOD=444-445 model for cloud
C     type ICLOUD.
C
C     Input variables
C	ICLOUD	INTEGER	Required cloud type
C	NVAR	INTEGER	Number of variable types in state vector
C       VARIDENT(MVAR,3) INTEGER Variables definition codes
C	VARPARAM(MVAR,MPARAM) REAL Associated parameter definition values
C	XN(MX)	REAL	State Vector
C	NPRO	INTEGER	Number of levels in vertical profile
C
C     Output variables
C	GETRADIUS REAL	radius of particles
C
C     Pat Irwin 25/5/16
C
C     **************************************************************
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      
      REAL XN(MX),X
      INTEGER NVAR,VARIDENT(MVAR,3),IVAR,NXTEMP,NP,NPVAR,NPRO
      INTEGER ICLOUD,J
      REAL VARPARAM(MVAR,MPARAM)

      NXTEMP=0
      X=-1.0
      DO 10 IVAR=1,NVAR
C       print*,'A',ICLOUD,IVAR,(VARIDENT(IVAR,J),J=1,3)
       IF(VARIDENT(IVAR,1).EQ.444.AND.VARIDENT(IVAR,2).EQ.ICLOUD)THEN
        X=EXP(XN(NXTEMP+1))
C        print*,'X = ',X
       ENDIF
       IF(VARIDENT(IVAR,1).EQ.445.AND.VARIDENT(IVAR,2).EQ.ICLOUD)THEN
        X=EXP(XN(NXTEMP+1))
C        print*,'X = ',X
       ENDIF
       IF(VARIDENT(IVAR,1).LE.100)THEN
        NP=NPVAR(VARIDENT(IVAR,3),NPRO,VARPARAM(IVAR,1))
C        print*,'D',(VARIDENT(IVAR,J),J=1,3),NP
       ELSE
        NP=1
        IF(VARIDENT(IVAR,1).EQ.888)NP=INT(VARPARAM(IVAR,1))
        IF(VARIDENT(IVAR,1).EQ.887)NP=INT(VARPARAM(IVAR,1))
        if(varident(ivar,1).eq.444)then
         if(varparam(ivar,2).gt.0.0)then
          np = 2+int(varparam(ivar,1))
         else
          np = 3
         endif
        endif
        IF(VARIDENT(IVAR,1).EQ.445)NP=3+(2*INT(VARPARAM(IVAR,1)))
        IF(VARIDENT(IVAR,1).EQ.222)NP=8
        IF(VARIDENT(IVAR,1).EQ.223)NP=9
        IF(VARIDENT(IVAR,1).EQ.224)NP=9
        IF(VARIDENT(IVAR,1).EQ.225)NP=11
        IF(VARIDENT(IVAR,1).EQ.226)NP=8
        IF(VARIDENT(IVAR,1).EQ.227)NP=7
       ENDIF

C       print*,'B',NP,XN(NXTEMP+1:NXTEMP+NP)
       NXTEMP=NXTEMP+NP

10    CONTINUE

      GETRADIUS=X

      RETURN

      END
