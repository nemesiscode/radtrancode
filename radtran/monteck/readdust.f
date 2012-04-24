      SUBROUTINE READDUST(IPFILE,H,NPRO,NCONT,DUST)
C     *************************************************************
C     Subroutine to read in a dust profile.
C
C     Input variables
C       IPFILE          CHARACTER*100 File name
C	H(MAXPRO)	REAL	Height grid of accompanying T/P profile
C	NPRO		INTEGER	Number of points in that profile
C
C     Output variables:
C	NCONT		INTEGER	Number of dust types
C	DUST(MAXPRO,MAXCON) REAL	Dust profiles
C
C     Pat Irwin		1/4/99
C
C     *************************************************************

      INTEGER NPRO,NCONT
      INCLUDE '../includes/arrdef.f'

      REAL DUST(MAXPRO,MAXCON),DUSTH(MAXPRO),H(MAXPRO)
      REAL TDUST(MAXPRO),XDUST,XH,SDUST(MAXPRO)
      CHARACTER*100 IPFILE,BUFFER

      CALL REMSP(IPFILE)
      CALL FILE(IPFILE,IPFILE,'prf')
      CALL LOCASE(IPFILE)

      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
1     FORMAT(A)
54    READ(1,1)BUFFER  
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)NN, NCONT
      DO 105 J = 1, NN
       READ(1,*)DUSTH(J),(DUST(J,I),I=1,NCONT)
C       PRINT*,DUSTH(J),(DUST(J,I),I=1,NCONT)
105   CONTINUE
      IF(J.LT.2)THEN
       WRITE(*,108)
108    FORMAT(' error reading profile')
       GOTO 119
      END IF
C     now sorting
      DO 114 K=1,NN 
      DO 114 I=1,NN-1   
      IF(ABS(DUSTH(I)-DUSTH(I+1)).LT.0.01)THEN
       WRITE(*,115)
115    FORMAT(' identical heights found')
       PRINT*,I,DUSTH(I),DUSTH(I+1)
       GOTO 119
      END IF
      IF(DUSTH(I).GT.DUSTH(I+1))THEN
       TD1=DUSTH(I+1)
       DUSTH(I+1)=DUSTH(I)
       DUSTH(I)=TD1
       DO 233 L=1,NCONT
        TD1=DUST(I+1,L)
        DUST(I+1,L)=DUST(I,L)
        DUST(I,L)=TD1
233    CONTINUE
      END IF
114   CONTINUE
C     now interpolating the input array to find values at profile
C     heights
       
      DO 244 K=1,NCONT 
        
       DO 243 I=1,NN
        TDUST(I)=DUST(I,K)
243    CONTINUE

       DO 203 I=1,NPRO
        SDUST(I)=0.0
203    CONTINUE
      
       DO 107 I=1,NPRO
          XH = H(I)
          CALL VERINT(DUSTH,TDUST,NN,XDUST,XH)
          DUST(I,K)=XDUST
107    CONTINUE
244   CONTINUE


109   CLOSE(UNIT=1) 
        
      WRITE(*,34) NCONT
34    FORMAT(' read in model with',I5,' aerosol types')
C
      RETURN   

119   CLOSE(UNIT=1)

      RETURN

      END

   
