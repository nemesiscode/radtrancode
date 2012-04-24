      SUBROUTINE WRLBLD
C     $Id: wrlbld.f,v 1.4 2011-06-23 09:11:21 irwin Exp $
C----------------------------------------------------------------------------
C_TITLE:  WRLBLD: to output the GENLBL arrays and Calculation record in
C         driver file format
C
C_ARGS:
C
C_KEYS:   PROG,ATMO,SPEC,VMS
C
C_DESCR:  
C
C_FILES:  UNIT 4: driving file output
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:    6feb87 SBC ORIGINAL VERSION
C         23jun92 SBC minor format mods for version 0
C         18nov92 SBC changed to use EMTEMP and BASET
C          1dec92 SBC changed  layer and path data to column format
C         23feb93 SBC fixed indeces on writing CONT
C                     and now only prompts for driver file if not set,
C                     for batch mode use.
C------------------------------------------------------------------------------
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/laycom.f'
C     local variables (used only in this subroutine)
      INTEGER I,J,K
C
      IF(NLAYER.EQ.0)THEN
        CALL WTEXT('no layers defined - nothing written')
        STOP
        END IF
      IF(NPATH.EQ.0)THEN
        CALL WTEXT('no paths defined - nothing written')
        STOP
        END IF
C     checking to see if there are any unused layers and removing if reqested
      I=1
101   CONTINUE
      DO 102 J=1,NPATH
      IF(IMOD(J).EQ.8)GOTO 102
      IF(IMOD(J).EQ.21)THEN
       IF(CLRLAY)THEN
        print*,'Must not clear layers for scattering limb calculation'
        print*,'Setting to NOCLRLAY'
       ENDIF
       CLRLAY=.FALSE.
      ENDIF
      DO 108 K=1,NLAYIN(J)
      IF(LAYINC(K,J).EQ.I)THEN
        I=I+1
        IF(I.LE.NLAYER)GOTO 101
        GOTO 107
        END IF
108   CONTINUE
102   CONTINUE
      IF(CLRLAY)THEN
C       removing layer I
        DO 103 J=I,NLAYER-1
        TEMP(J)=TEMP(J+1)
        PRESS(J)=PRESS(J+1)
        HFP(J)=HFP(J+1)
        DOP(J)=DOP(J+1)
        BASEP(J)=BASEP(J+1)
        BASEH(J)=BASEH(J+1)
        DELH(J)=DELH(J+1)
        BASET(J)=BASET(J+1)
        TOTAM(J)=TOTAM(J+1)
        DO 104 K=1,NGAS
        AMOUNT(J,K)=AMOUNT(J+1,K)
        PP(J,K)=PP(J+1,K)
104     CONTINUE
        DO 105 K=1,NCONT
        CONT(K,J)=CONT(K,J+1)
105     CONTINUE
103     CONTINUE
        DO 109 J=1,NPATH
C       don't need to remove layer reference for combined paths
        IF(IMOD(J).EQ.7.OR.IMOD(J).EQ.8)GOTO 109
        DO 106 K=1,NLAYIN(J)
        IF(LAYINC(K,J).GT.I)LAYINC(K,J)=LAYINC(K,J)-1
106     CONTINUE
109     CONTINUE
        NLAYER=NLAYER-1
        GOTO 101
        END IF
107     CONTINUE

C
      CALL REMSP(DRFILE)
      IF(DRFILE(1:1).EQ.' ')THEN
        CALL PROMPT('name of driving file?')
        READ(*,503)OPFILE
503     FORMAT(A)
        CALL FILE(OPFILE,DRFILE,'drv')
       END IF
      OPEN(UNIT=4,FILE=DRFILE,STATUS='UNKNOWN')
      CALL FILE(DRFILE,OPFILE,'out')
C
C     output section
      WRITE(4,516)DRFILE
516   FORMAT(' original name of this file: ',1A30)
      IF(ICONV.GE.0.AND.ICONV.LT.10)THEN
       WRITE(4,507)VMIN,DELV,NPOINT,FWHM
507    FORMAT(1X,F10.3,F10.5,I10,F8.3,4X,' :Vmin dV Npoints FWHM-LBL')
       WRITE(4,508)WING,VREL
508    FORMAT(1X,2F8.3,6X,' : wing continuum limit and overall limit')
      ELSE IF(ICONV.GE.10.AND.ICONV.LT.20)THEN
       WRITE(4,607)VMIN,DELV,NPOINT,FWHM
607    FORMAT(1X,F10.3,F10.5,I10,F8.3,4X,' :Vmin dV Npoints FWHM-BAND')
       WRITE(4,608)WING,VREL
608    FORMAT(1X,2F8.3,6X,' : Additional codes PAR1 and PAR2')
      ELSE IF(ICONV.GE.20.AND.ICONV.LT.30)THEN
       WRITE(4,707)VMIN,DELV,NPOINT,FWHM
707    FORMAT(1X,F10.3,F10.5,I10,F8.3,4X,' :Vmin dV Npts FWHM-CORRK')
       WRITE(4,708)WING,VREL
708    FORMAT(1X,2F8.3,6X,' : Additional codes PAR1 and PAR2')
      ELSE
       WRITE(*,*)'WRLBLD: ICONV out of range'
       STOP
      END IF
      WRITE(4,560)LINKEY
560   FORMAT(1X,A)
      WRITE(4,509)ICONV,FLAGH2P,NCONT,FLAGC
509   FORMAT(1X,I4,I4,I6,I4,12X,
     1 ' : spectral model code, FLAGH2P, NCONT, FLAGC')
      IF(NCONT.GT.0)WRITE(4,530)XFILE
530   FORMAT(1X,A30,' : Dust x-section file')
      WRITE(4,501)NLAYER,NPATH,NGAS
501   FORMAT(1X,3I4,10X,' : number of layers, paths and gases')
      DO 504 I=1,NGAS
      WRITE(4,510)IDGAS(I),I
      WRITE(4,511)ISOGAS(I),IPROC(I)
504   CONTINUE
510   FORMAT(1X,I4,35X,' : identifier for gas',I2)
511   FORMAT(1X,I6,2X,I2,29X,' : isotope ID and process parameter')
C     outputting details of each layer
      WRITE(4,551)
551   FORMAT('format of layer data')
      WRITE(4,552)
552   FORMAT(
     1' layer baseH  delH   baseP      baseT   ',
     1'totam       pressure    temp   doppler')
      WRITE(4,553)
553   FORMAT('        absorber amounts and partial pressures')
      WRITE(4,554)
554   FORMAT('        continuum points if any')
      DO 505 I=1,NLAYER
      WRITE(4,512)I,BASEH(I),DELH(I),BASEP(I),BASET(I),
     1TOTAM(I),PRESS(I),TEMP(I),DOP(I)
512   FORMAT(I3,1X,F7.1,1X,F7.1,1X,E11.5,1X,F8.3,1X,E11.5,
     11X,E11.5,1X,F8.3,1X,F7.4)
      WRITE(4,513)(AMOUNT(I,J),PP(I,J),J=1,NGAS)
513   FORMAT(8X,6(1X,E11.5))
514   FORMAT(8X,1X,E11.5,10(1X,I3))
C      IF(NCONT.GT.1)WRITE(4,513)(CONT(J,I),J=1,NCONT)
      IF(NCONT.GT.0)WRITE(4,513)(CONT(J,I),J=1,NCONT)
      IF(FLAGH2P.EQ.1)WRITE(4,513)HFP(I)
      IF(FLAGC.EQ.1)WRITE(4,514)HFC(I),(IFC(J,I),J=1,NCONT)
505   CONTINUE
      DO 520 I=1,NPATH
      WRITE(4,521)NLAYIN(I),IMOD(I),ERRLIM(I),I
521   FORMAT(1X,2I4,E12.5,12X,': Nlayers, model & error limit, path',I3)
      DO 522 J=1,NLAYIN(I)
      WRITE(4,523)J,LAYINC(J,I),EMTEMP(J,I),SCALE(J,I)
523   FORMAT(1X,I3,1X,I4,1X,F8.3,1X,E12.5,
     1'  :     layer or path, emission temp, scale')
522   CONTINUE
520   CONTINUE
      WRITE(4,525)NFILT
525   FORMAT(1X,I4,35X,' : number of filter profile points')
      DO 526 I=1,NFILT
      WRITE(4,527)FILTER(I),VFILT(I),I
527   FORMAT(1X,E12.5,F10.3,17X,' : filter profile point',I4)
526   CONTINUE
      WRITE(4,540)OPFILE
540   FORMAT(1X,A)
      WRITE(4,541)NCALC
541   FORMAT(1X,I3,37X,':number of calculations')
      DO 542 I=1,NCALC
      WRITE(4,543)ITYPE(I),NINTP(I),NREALP(I),NCHP(I),I
543   FORMAT(1X,4I4,24X,':type and # of parameters for calc',I3)
      DO 544 J=1,NINTP(I)
      WRITE(4,*)ICALD(J,I)
544   CONTINUE
      DO 545 J=1,NREALP(I)
      WRITE(4,*)RCALD(J,I)
545   CONTINUE
      DO 546 J=1,NCHP(I)
      WRITE(4,547)CCALD(J,I)
547   FORMAT(1X,1A30)
546   CONTINUE
542   CONTINUE

      CLOSE(4)

      RETURN
      END
C-----------------------------------------------------------------------------
