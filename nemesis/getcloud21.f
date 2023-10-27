      SUBROUTINE GETCLOUD21(ICLOUD,NVAR,VARIDENT,VARPARAM,XN,NPRO,
     1 XDEEP,HKNEE,XFSH,REFRADIUS)
C     **************************************************************
C     Function to sift through XN and VARIDENT FILES to find the mean 
C     radius of particles defined though IMOD=444-446 model for cloud
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
C	XDEEP REAL	Optical Depth
C	HKNEE REAL	Base cloud altitude
C	XFSH  REAL	Reference fractional scale height.
C
C     Pat Irwin 26/5/16
C
C     **************************************************************
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      
      REAL XN(MX),X
      INTEGER NVAR,VARIDENT(MVAR,3),IVAR,NXTEMP,NP,NPVAR,NPRO
      INTEGER ICLOUD,J
      REAL VARPARAM(MVAR,MPARAM),XDEEP,HKNEE,XFSH,REFRADIUS
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      NXTEMP=0
      XDEEP=-1.0
      HKNEE=-1.0
      XFSH=-1.0

      DO 10 IVAR=1,NVAR
C       if(idiag.gt.0)print*,'A',ICLOUD,IVAR,(VARIDENT(IVAR,J),J=1,3)
       IF(VARIDENT(IVAR,3).EQ.21.AND.-VARIDENT(IVAR,1).EQ.ICLOUD)THEN
        XDEEP=EXP(XN(NXTEMP+1))
        HKNEE=XN(NXTEMP+2)
        XFSH=VARPARAM(IVAR,1)
        REFRADIUS=VARPARAM(IVAR,2)
       ENDIF
       IF(VARIDENT(IVAR,1).LE.100)THEN
        NP=NPVAR(VARIDENT(IVAR,3),NPRO,VARPARAM(IVAR,1))
C        if(idiag.gt.0)print*,'D',(VARIDENT(IVAR,J),J=1,3),NP
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
        if(varident(ivar,1).eq.446)then
         if(varparam(ivar,2).gt.0.0)then
          np = 3+2*int(varparam(ivar,1))
         else
          np = 5
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

C       if(idiag.gt.0)print*,'B',NP,XN(NXTEMP+1:NXTEMP+NP)
       NXTEMP=NXTEMP+NP

10    CONTINUE

      RETURN

      END
