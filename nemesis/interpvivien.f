      SUBROUTINE INTERPVIVIEN(XLAT,XLON,NPRO,NVMR,P,H,T,VMR)

      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INTEGER NPRESS,NLAT,NLON,NGV
      PARAMETER (NPRESS=53,NLON=64,NLAT=32,NGV=6)
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL VP(NPRESS),VLON(NLON),VLAT(NLAT)
      REAL VT(NLON,NLAT,NPRESS),VVMR(NGV,NLON,NLAT,NPRESS)
      REAL XLAT,XLON,LP1
      REAL VY1(NPRESS),VY2(NPRESS),VY3(NPRESS),VY4(NPRESS)
      REAL Y1,Y2,Y3,Y4
      INTEGER NPRESS1,I,J,K,L,JLAT,JLON1,JLON2,IVMR
      INTEGER NLAT1,NLON1,NPRO,NVMR
      REAL FLAT,FLON
      CHARACTER*100 ANAME
      COMMON /VSTORE/ NLON1,NLAT1,NPRESS1,VP,VT,VVMR


      print*,'Calling interpvivien'
      print*,'NLON1 = ',NLON1
C      print*,'XLAT,XLON = ',XLAT,XLON
         
      IF(NLON1.NE.64)THEN
       ANAME='process_vivien.txt'
       OPEN(12,FILE=ANAME,STATUS='OLD')
        READ(12,*)NLON1,NLAT1
        IF(NLON.NE.NLON1)THEN
         PRINT*,'Error in interpvivien.f'
         PRINT*,'NLON1 <> NLON',NLON1,NLON
         STOP
        ENDIF
        IF(NLAT.NE.NLAT1)THEN
         PRINT*,'Error in interpvivien.f'
         PRINT*,'NLAT1 <> NLAT',NLAT1,NLAT
         STOP
        ENDIF
        READ(12,*)(VLON(I),I=1,NLON)
        READ(12,*)(VLAT(I),I=1,NLAT)
        READ(12,*)NPRESS1
        IF(NPRESS.NE.NPRESS1)THEN
         PRINT*,'Error in interpvivien.f'
         PRINT*,'NPRESS1 <> NPRESS',NPRESS1,NPRESS
         STOP
        ENDIF

        READ(12,*)(VP(I),I=1,NPRESS)
  
        DO I=1,NLON
         DO J=1,NLAT
          DO K=1,NPRESS
           READ(12,*)VT(I,J,K),(VVMR(L,I,J,K),L=1,6) 
          ENDDO
         ENDDO
        ENDDO
  
       CLOSE(12)

C      Convert P from bar to atm and convert to log
       DO I=1,NPRESS
        VP(I)=ALOG(VP(I)/1.013)
       ENDDO

C      Vivien has 0 longitude as sub-stellar point. He we have 180
C      Hence need to add 180 to all longitudes
       DO I=1,NLON
        VLON(I)=VLON(I)+180.0
       ENDDO

      ENDIF



C     Find closest point in stored array

      JLAT=-1
      DO I=1,NLAT-1
       IF(XLAT.GE.VLAT(I).AND.XLAT.LE.VLAT(I+1))THEN
        JLAT=I
        FLAT = (XLAT-VLAT(I))/(VLAT(I+1)-VLAT(I))
       ENDIF
      ENDDO
      IF(JLAT.LT.0)THEN
        IF(XLAT.LT.VLAT(1))THEN
         JLAT=1
         FLAT=0.
        ENDIF
        IF(XLAT.GE.VLAT(NLAT))THEN
         JLAT=NLAT-1
         FLAT=1.0
        ENDIF
      ENDIF

      IF(JLAT.LT.0.0)THEN
        PRINT*,'Cannot locate JLAT'
        PRINT*,XLAT,VLAT(1),VLAT(NLAT)
        STOP
      ENDIF

      JLON1=-1
      JLON2=-1
      DO I=1,NLON-1
       IF(XLON.GE.VLON(I).AND.XLON.LE.VLON(I+1))THEN
        JLON1=I
        JLON2=I+1
        FLON = (XLON-VLON(I))/(VLON(I+1)-VLON(I))
       ENDIF
      ENDDO

      IF(JLON1.LT.0)THEN
       IF(XLON.LT.VLON(1))THEN
C        Xlon must be in range 0. to VLON(1)
         JLON1=NLON
         JLON2=1
         FLON=(XLON+360.0-VLON(NLON))/(VLON(1)+360.0-VLON(NLON))
       ENDIF
       IF(XLON.GE.VLON(NLON))THEN
C        Xlon must be in range VLON(NLON) to 360.
         JLON1=NLON
         JLON2=1
         FLON=(XLON-VLON(NLON))/(VLON(1)+360.0-VLON(NLON))
       ENDIF
      ENDIF

      IF(JLON1.LT.0.0)THEN
        PRINT*,'Cannot locate JLON1'
        PRINT*,XLON,VLON(1),VLON(NLON)
        STOP
      ENDIF

C      PRINT*,'JLAT,FLAT = ',JLAT,FLAT
C      PRINT*,'VLAT(JLAT),VLAT(JLAT+1),XLAT',VLAT(JLAT),VLAT(JLAT+1),XLAT
C      PRINT*,'JLON1,JLON2,FLON',JLON1,JLON2,FLON
C      PRINT*,'VLON(JLON1),VLON(JLON2),XLON',VLON(JLON1),VLON(JLON2),XLON

      DO I=1,NPRESS
       VY1(I)=VT(JLON1,JLAT,I)
       VY2(I)=VT(JLON2,JLAT,I)
       VY3(I)=VT(JLON2,JLAT+1,I)
       VY4(I)=VT(JLON1,JLAT+1,I)
      ENDDO

      DO I=1,NPRO
       LP1=ALOG(P(I))
       CALL VERINT(VP,VY1,NPRESS,Y1,LP1)
       CALL VERINT(VP,VY2,NPRESS,Y2,LP1)
       CALL VERINT(VP,VY3,NPRESS,Y3,LP1)
       CALL VERINT(VP,VY4,NPRESS,Y4,LP1)
       T(I) = (1.0-FLON)*(1.0-FLAT)*Y1 + FLON*(1.0-FLAT)*Y2+
     1 FLON*FLAT*Y3 + (1.0-FLON)*FLAT*Y4
      ENDDO

      DO IVMR=1,NVMR

       DO I=1,NPRESS
        VY1(I)=VVMR(IVMR,JLON1,JLAT,I)
        VY2(I)=VVMR(IVMR,JLON2,JLAT,I)
        VY3(I)=VVMR(IVMR,JLON2,JLAT+1,I)
        VY4(I)=VVMR(IVMR,JLON1,JLAT+1,I)
       ENDDO

       DO I=1,NPRO
        LP1=ALOG(P(I))
        CALL VERINT(VP,VY1,NPRESS,Y1,LP1)
        CALL VERINT(VP,VY2,NPRESS,Y2,LP1)
        CALL VERINT(VP,VY3,NPRESS,Y3,LP1)
        CALL VERINT(VP,VY4,NPRESS,Y4,LP1)
        VMR(I,IVMR) = (1.0-FLON)*(1.0-FLAT)*Y1 + FLON*(1.0-FLAT)*Y2+
     1 FLON*FLAT*Y3 + (1.0-FLON)*FLAT*Y4

       ENDDO

      ENDDO

      RETURN

      END
