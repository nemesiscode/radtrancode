      SUBROUTINE DISPAT
C     $Id: display.f,v 1.2 2011-06-17 15:40:26 irwin Exp $
C----------------------------------------------------------------------------
C_TITLE:  DISPAT:
C
C_ARGS:
C
C_KEYS:   SUBR,ATMO,SPEC,VMS
C
C_DESCR:  Outputs path data to screen  
C
C_FILES:
C
C_CALLS:  Called by Pl_spec
C
C_BUGS:
C
C_HIST:   10feb87 SBC  ORIGINAL VERSION
C	  26apr12 PGJI Added comments and updated for new radtrancode suite.
C
C------------------------------------------------------------------------------
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'

C-----------------------------------------------------------------------------
C     local variables (used only in this program)
      INTEGER I,J
C-----------------------------------------------------------------------------
      DO 701 I=1,NPATH
      WRITE(*,702)I,IMOD(I),NLAYIN(I)
702   FORMAT(' path:',I2,'  model type:',I2,',  ',I2,'layers included')
      WRITE(*,703)(LAYINC(J,I),J=1,NLAYIN(I))
703   FORMAT(1X,25I3)
      WRITE(*,704)
704   FORMAT(' ---------------------------------------------------')
701   CONTINUE
      RETURN
      END
C
C------------------------------------------------------------------------------
      SUBROUTINE DISLAY
C----------------------------------------------------------------------------
C_TITLE:  DISLAY:
C
C_ARGS:
C
C_KEYS:   SUBR,ATMO,SPEC,VMS
C
C_DESCR:  
C
C_FILES:
C
C_CALLS:
C
C_BUGS:
C
C_HIST:   10feb87 SBC ORIGINAL VERSION
C------------------------------------------------------------------------------
C     note that include is not F77 and is included only during development
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C-----------------------------------------------------------------------------
C     local variables (used only in this program)
      INTEGER I,J
C-----------------------------------------------------------------------------
      WRITE(*,806)NGAS
806   FORMAT(' there are',I3,' gases')
      DO 807 I=1,NGAS
      WRITE(*,808)I,IDGAS(I),ISOGAS(I)
808   FORMAT(' gas:',I2,'  identifier:',I3,'   isotope id:',I5)
807   CONTINUE
      WRITE(*,805)
      DO 801 I=1,NLAYER
      WRITE(*,802)I,PRESS(I),TEMP(I)
802   FORMAT(' layer:',I2,'  pressure:',E12.5,'atm  temp:',F7.2,'K')
      WRITE(*,804)
804   FORMAT(' gas amounts')
      WRITE(*,803)(AMOUNT(I,J),J=1,NGAS)
803   FORMAT(1X,6E12.5)
      WRITE(*,805)
805   FORMAT(' --------------------------------------------------')
801   CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------------
      SUBROUTINE DISCAL
C-----------------------------------------------------------------------------
C_TITLE:  DISCAL
C
C_ARGS:
C
C_KEYS:   SUBR,ATMO,SPEC,VMS
C
C_DESCR:  
C
C_FILES:
C
C_CALLS:
C
C_BUGS:
C
C_HIST:   10feb87 SBC ORIGINAL VERSION
C------------------------------------------------------------------------------
C     note that include is not F77 and is included only during development
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C-----------------------------------------------------------------------------
C     local variables (used only in this program)
      INTEGER I,J
      LOGICAL BIT
      CHARACTER*40 TEXT
C-----------------------------------------------------------------------------
      WRITE(*,901)NCALC
901   FORMAT(' there are',I3,' calculations')
      DO 902 I=1,NCALC
      WRITE(*,*)' '
      WRITE(*,903)
903   FORMAT(' --------------------------------------------------')
      WRITE(*,866)I
866   FORMAT(' Calculation : ',I2)
      CALL DISCAL1(I)
902   CONTINUE
      WRITE(*,903)
      WRITE(*,*)' '
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE DISCAL1(I)
C-----------------------------------------------------------------------------
C_TITLE:  DISCAL1
C
C_ARGS:
C
C_KEYS:   SUBR,ATMO,SPEC,VMS
C
C_DESCR:  
C
C_FILES:
C
C_CALLS:
C
C_BUGS:
C
C_HIST:   10feb87 SBC ORIGINAL VERSION
C	  21oct94 PGJI modified from DISCAL to display single calculation
C------------------------------------------------------------------------------
C     note that include is not F77 and is included only during development
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C-----------------------------------------------------------------------------
C     local variables (used only in this program)
      INTEGER I,J,FSTCEL
      LOGICAL BIT
      CHARACTER*40 TEXT
C-----------------------------------------------------------------------------

      IF(ITYPE(I).EQ.160.OR.ITYPE(I).EQ.161)THEN
       WRITE(*,*)'Combined atmosphere and cell calculations'
       IF(ITYPE(I).EQ.161)WRITE(*,*)'Reflecting layer transmissions'
       NPATH1=ICALD(1,I)
       NPATH2=ICALD(2,I)
       IC=ICALD(3,I)
       IA=ICALD(4,I)
       IF(NPATH2.EQ.NPATH1)THEN
        WRITE(*,930)NPATH2
       ELSE
        WRITE(*,931)NPATH1,NPATH2
       END IF
930    FORMAT(' Combined output path is : ',I3)
931    FORMAT(' Output paths run from: ',I3,' to ',I3)
       WRITE(*,932)IC
       WRITE(*,933)IA
932    FORMAT(' Cell calculation number is : ',I3)
933    FORMAT(' Atmosphere calculation number is : ',I3)
      ELSE IF(ITYPE(I).GE.128.AND.ITYPE(I).LE.131)THEN
       WRITE(*,*)'Cell calculation'
       IF(ITYPE(I).EQ.128)THEN
        FSTCEL=ICALD(1,I)
        LSTCEL=ICALD(2,I)
        NSNGL=1+FSTCEL-LSTCEL
        WRITE(*,*)'Single cell calculations'
        WRITE(*,300)NSNGL
300     FORMAT(' Number of cells = ',I3)
        DO 311 J=1,NSNGL
         LAY1 = LAYINC(1,I)
         IP1 = FSTCEL + (J-1)
         WRITE(*,301)IP1,LAY1
301      FORMAT(' Path: ',I3,' Layer: ',I3)
311     CONTINUE
       ELSE IF(ITYPE(I).EQ.129)THEN
        WRITE(*,*)'PMR - 2 pressure approximation'
        WRITE(*,312)ICALD(1,I)
312     FORMAT(' Sideband transmission. Path : ',I3)
        WRITE(*,313)ICALD(2,I)
313     FORMAT(' Wideband transmission. Path : ',I3)
       ELSE IF(ITYPE(I).EQ.130)THEN
        WRITE(*,*)'PMR File '
        NPHAS=NREALP(I)
        FSTCEL=ICALD(1,I)
        LSTCEL=ICALD(2,I)
        WRITE(*,314)NPHAS
314     FORMAT(' Number of points around cycle : ',I3)
        DO 315 J=1,NPHAS
         ICEL=FSTCEL + (J-1)
         WRITE(*,316)J,ICEL,RCALD(J,I)
316      FORMAT(' J = ',I3,' Path = ',I3,' Phase = ',F8.2)
315     CONTINUE
       ELSE IF(ITYPE(I).EQ.131)THEN
        WRITE(*,*)'SCR'
        WRITE(*,312)ICALD(1,I)
        WRITE(*,313)ICALD(2,I)
       END IF

      ELSE

       WRITE(*,*)'Atmosphere calculation'
       IF(ITYPE(I).EQ.256)THEN
        WRITE(*,*)'Scattering calculation'
       ELSE
        IF(BIT(4,ITYPE(I)))THEN
         WRITE(*,*)'Single atmospheric path'
        END IF
        IF(BIT(5,ITYPE(I)))THEN
         WRITE(*,*)'Reflecting layer atmosphere transmission'
        ELSE IF(BIT(6,ITYPE(I)))THEN
         WRITE(*,*)'Limb View'
        ELSE
         WRITE(*,*)'Nadir View'
        END IF
        IF(BIT(3,ITYPE(I)))THEN
         WRITE(*,*)'Curtis Godson approximation has been used'
        END IF
        IF(BIT(2,ITYPE(I)))THEN
         WRITE(*,*)'Weighting function calculation'
         GOTO 222
        END IF
        IF(BIT(1,ITYPE(I)))THEN
         WRITE(*,*)'Thermal emission calculation'
         GOTO 222
        END IF
        IF(BIT(0,ITYPE(I)))THEN
         WRITE(*,*)'Absorption calculation'
        ELSE
         WRITE(*,*)'Transmission calculation'
        END IF
       END IF

222    CONTINUE
       IF(.NOT.BIT(4,ITYPE(I)))THEN
        IF(.NOT.BIT(5,ITYPE(I)))THEN
        WRITE(*,223)RCALD(1,I)
223     FORMAT(' Viewing angle = ',F8.2)
        WRITE(*,224)RCALD(2,I)
224     FORMAT(' Base height of lowest layer (km) = ',F8.2)
        NPATH1=ICALD(1,I)
        NPATH2=ICALD(2,I)
       ELSE
        WRITE(*,243)RCALD(1,I),RCALD(2,I)
243     FORMAT(' Incident and reflected angles = ',F8.2,' ',F8.2)
        WRITE(*,224)RCALD(3,I)
        NPATH1=ICALD(1,I)
        NPATH2=ICALD(2,I)
       END IF
       IF(NPATH2.EQ.NPATH1)THEN
        WRITE(*,225)NPATH1
       ELSE
        WRITE(*,226)NPATH1,NPATH2
       END IF
225    FORMAT(' 1 atmospheric path  = ',I3)
226    FORMAT(' Calculation involves paths: ',I3,' to ',I3)
       END IF
      END IF

      IF(ITYPE(I).EQ.200)THEN
       WRITE(*,*)'Reflecting Layer Calculation'
       WRITE(*,*)'Pure atmospheric paths'
       KE=ICALD(1,I)
       KR=ICALD(2,I)
       WRITE(*,366)KE
       WRITE(*,356)KR
366    FORMAT('Number of atmosphere calculations = ',I3)
356    FORMAT('Number of reflected atmosphere calcs = ',I3)
       WRITE(*,*)'Atmospheric calculations are:'
       DO K=1,KE
        WRITE(*,*)ICALD(2+K,I)
       END DO
       WRITE(*,*)'Reflection Atmosphere calculations are:'
       DO K=1,KR
        WRITE(*,*)ICALD(2+K+KE,I)
       END DO
      END IF
      IF(ITYPE(I).EQ.201)THEN
       WRITE(*,*)'Reflecting Layer Calculation'
       WRITE(*,*)'Combined cell/atmospheric paths'
       KE=ICALD(1,I)
       KR=ICALD(2,I)
       WRITE(*,366)KE
       WRITE(*,356)KR
       WRITE(*,*)'Cell/Atmospheric calculations are:'
       DO K=1,KE
        WRITE(*,*)ICALD(2+K,I)
       END DO
       WRITE(*,*)'Cell/Reflection Atmosphere calculations are:'
       DO K=1,KR
        WRITE(*,*)ICALD(2+K+KE,I)
       END DO
      END IF


      WRITE(*,*)' '
      RETURN
      END
C-----------------------------------------------------------------------------
