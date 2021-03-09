      SUBROUTINE RFCMOD(TEXT)
C     $Id: rfcmod.f,v 1.1 2007-06-28 15:02:35 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  RFPMOD: reads in fractional cloud cover profile for path.f
C
C_KEYS:   RADTRAN,SUBR
C
C_DESCR:  
C
C_ARGS:   
C
C_FILES : unit 1 - atmospheric model, dust and cell input files
C         unit 2 - the path file [.pat]
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   03Dec93 PGJI Original version
C
C_END:
C--------------------------------------------------------------
      CHARACTER TEXT*(*)
C--------------------------------------------------------------
C     Variables to hold calculated layers and the details of each paths and
C     calculation requested.
C     pathcom holds the variables used bu the genlbl software too
C     laycom holds variables used only by the path software
C     parameters are passed between routines mostly using common blocks
C     because of the extensive use of large arrays
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/laycom.f'
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in code
      CHARACTER*100 IPFILE,BUFFER
      INTEGER I,J,K
      REAL TD1,TMP(500)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet
C--------------------------------------------------------------
C
C     reading in in vertical profiles produced by profile.for
C     path file format is:
C        model filename
C
      READ(TEXT,1)IPFILE
1     FORMAT(A)
      CALL REMSP(IPFILE)
      CALL LOCASE(IPFILE)
      CALL FILE(IPFILE,IPFILE,'prf')
      WRITE(*,*)' RFCMOD.f :: reading fcover-model: ',ipfile
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
54    READ(1,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)NN,NCONT1
      DO 105 J = 1, NN     
       READ(1,*)FCLOUDH(J),FCLOUD(J),(ICLOUD(I,J),I=1,NCONT1)
105   CONTINUE
      CLOSE(1)

      IF(J.LT.2)THEN
       WRITE(*,108)
108    FORMAT(' error reading FCLOUD profile')
       CLOSE(1)
       STOP
      END IF

C     now sorting
      DO 114 K=1,NN
      DO 112 I=1,NN-1
      IF(ABS(FCLOUDH(I)-FCLOUDH(I+1)).LT.0.01)THEN
       WRITE(*,115)
115    FORMAT('rfcmod: identical heights found')
       print*,K,I,FCLOUDH(I),FCLOUDH(I+1)
       STOP
      END IF
      IF(FCLOUDH(I).GT.FCLOUDH(I+1))THEN
       if(idiag.gt.0)then
        print*,'rfpmod: reordering'
        print*,I,FCLOUDH(I),FCLOUDH(I+1)
       endif
       TD1=FCLOUDH(I+1)
       FPH2H(I+1)=FCLOUDH(I)
       FCLOUDH(I)=TD1
       TD1=FCLOUD(I+1)
       FCLOUD(I+1)=FCLOUD(I)
       FCLOUD(I)=TD1
      END IF
112   CONTINUE
114   CONTINUE
C     now interpolating the input array to find values at profile
C     heights

      DO 107 I=1,NPRO
          CALL VERINT(FCLOUDH,FCLOUD,NN,FCLOUDI(I),H(I))
107   CONTINUE

      DO 109 J= 1,NCONT1
            DO K=1,NN
             TMP(K)=1.0*ICLOUD(J,K)
            ENDDO
            DO I=1,NPRO
             CALL VERINT(FCLOUDH,TMP,NN,XOUT,H(I))
             ICLOUDI(J,I)=INT(XOUT+0.5)
            ENDDO
109   CONTINUE

      JFC = 1

      RETURN

      END
