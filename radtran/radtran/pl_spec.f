           PROGRAM PLSPEC
C     $Id: pl_spec.f,v 1.8 2011-01-20 14:34:51 teanby Exp $
C-----------------------------------------------------------------------------
C_TITLE:  PLTSPEC: to plot output of GENLBL calculations
C
C_ARGS:   none
C
C_KEYS:   ATMO,SPEC,VMS,PROG
C
C_DESCR:  plots any one of the spectra in a file produced by LBL
C
C_CALLS:
C
C_BUGS:
C
C_HIST:   5feb87  SBC ORIGINAL VERSION
C         8MAY87  SBC Modified for GKS version of NCAR graphics
C         27aug87 SBC Modified to use SIMPLEPLOT
C         23jun92 SBC hacked out of plot_spec as simple spectrum only routine
C         21oct94 PGJI Major revisions:
C			1)DISCAL made more informative
C			2)Calculations section completed for everything except
C			PMR Files and SCRs
C			3)Weighting function plotting added
C-----------------------------------------------------------------------------
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C-----------------------------------------------------------------------------
      INTEGER MAXPT
      PARAMETER(MAXPT=150000)
      INTEGER IREC,IREC0,IPATH,ICALC,I,J,NPLOT,IPLOT,NINTER,NCHX,NCHY
      REAL VMAX,YMIN,YMAX,XMIN,XMAX,VPLOT,OUTPUT,DX,DY,OLDOUT
      REAL Y(MAXPT),Y1(100),BRIGHT,PLOG10,ERROR(MAXPT)
      double precision X(MAXPT),ddelv,dvmin
      character*20 tmp
      REAL ERROUT
      CHARACTER*100 XLABEL,YLABEL,FILE2
      CHARACTER Q
      CHARACTER*100 IPFILE,XFORM,YFORM
      INTEGER ISPACE,LSTCEL,FSTCEL,CELLPATH
      LOGICAL ABSORB,FIRST,LOG,BIT,WEIGHTF,ASCOUT,EMISS,MICRON
      LOGICAL BRTEMP,COMPLE
      INTEGER RPATH(10),APATH(10)
      REAL ALBEDO(10),PI,THETA0,SOLZEN(10)
      PARAMETER(PI=3.1415927,THETA0=0.5*9.30483E-03)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      idiag=1

      ABSORB=.FALSE.
      BRTEMP=.FALSE.
      
C     *** start of debug lines ***
C      PRINT*,'=== TEST Pl_spec ==='
C     *** end of debug lines ***

      CALL PROMPT('driving file name?')
      READ(*,503)OPFILE
503   FORMAT(1A30)
      CALL FILE(OPFILE,IPFILE,'drv')
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      CALL RDLBLD
C
C     check if want to look at averaged file
      CALL PROMPT('averaged?')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        CALL FILE(OPFILE,FILE2,'ave')
       ELSE
        FILE2=OPFILE
      END IF
      I=2*ISYS()
      PRINT*,'Record length = ',I
      OPEN(UNIT=2,FILE=FILE2,STATUS='OLD',ACCESS='DIRECT',RECL=I)
      READ(2,REC=1)IREC0
      READ(2,REC=2)NPOINT
      READ(2,REC=3)VMIN
      READ(2,REC=4)DELV
      READ(2,REC=5)NPATH

C-----------------------------------------------------
c	some fudge to fix the wavnumber scale for very small spacings! 	
C-----------------------------------------------------
      open(66,file='.tmp_pl_spec',status='unknown')
      write(66,*) vmin,delv
      close(66)
	open(66,file='.tmp_pl_spec',status='old')
      read(66,*) dvmin,ddelv
      close(66)
	print*,dvmin,ddelv
C-----------------------------------------------------
      
      PRINT*,'header: IREC0,NPOINT,VMIN,DELV,NPATH',
     1  IREC0,NPOINT,VMIN,DELV,NPATH

      IF(NPOINT.GT.MAXPT)THEN
        WRITE(*,*)'NPOINT',NPOINT,'MAXPT',MAXPT
        WRITE(*,301)
301     FORMAT('Pl_spec: NPOINT > MAXPT. Abort')
        STOP
      END IF

C      Debugging lines
C      NOUT = NPATH*NPOINT
C      IREC=IREC0
C      DO I=1,NOUT
C       READ(2,REC=IREC)Y(I),ERROR(I)
C       print*,IREC,Y(I),ERROR(I)
C       IREC=IREC+1
C      END DO

C
C-----------------------------------------------------
10    CONTINUE
      WRITE(*,11)
      WRITE(*,12)
      WRITE(*,13)
      WRITE(*,14)
      WRITE(*,15)
      WRITE(*,16)
11    FORMAT(' A - display calculations')
12    FORMAT(' B - display paths')
13    FORMAT(' C - display layers')
14    FORMAT(' D - plot a spectrum for a path')
15    FORMAT(' E - plot results of a calculation')
16    FORMAT(' Q - quit')
      CALL PROMPT('option')
      READ(*,21)Q
21    FORMAT(1A1)
      CALL UPCASE(Q)
      IF(Q.EQ.'A')THEN
        CALL DISCAL
       ELSE IF(Q.EQ.'B')THEN
        CALL DISPAT
       ELSE IF(Q.EQ.'C')THEN
        CALL DISLAY
       ELSE IF(Q.EQ.'D')THEN
          WRITE(*,*)'*** This option may not work correctly ***'
	  WRITE(*,*)'******** in some circumstances **********'
          IF(NPATH.EQ.1)THEN
             IPATH=1
             GOTO 131
          END IF
102     WRITE(*,101)NPATH
101     FORMAT(' there are ',I3,' paths')
        CALL PROMPT('plot which one?')
        READ(*,*)IPATH
        IF(IPATH.LT.1.OR.IPATH.GT.NPATH)GOTO 102
        GOTO 131
       ELSE IF(Q.EQ.'E')THEN
          IF(NCALC.EQ.1)THEN
             ICALC=1
             GOTO 203
          END IF
205     WRITE(*,207)NCALC
207     FORMAT('there are ',I3,' calculations')
        CALL PROMPT('plot which one?')
        READ(*,*)ICALC
        IF(ICALC.LT.1.OR.ICALC.GT.NCALC)GOTO 205
203     CONTINUE
        WRITE(*,*)' '
        WRITE(*,*)'------------------------------------------'
C     *** start of debug lines ***
C      PRINT*,'ITYPE(ICALC) = ',ITYPE(ICALC)
C     *** end of debug lines ***
        CALL DISCAL1(ICALC)
        WRITE(*,*)'------------------------------------------'
        WRITE(*,*)' '
        IF(ITYPE(ICALC).EQ.160.OR.ITYPE(ICALC).EQ.161)THEN
           WRITE(*,*)'Combined Atmosphere and cell Calculation'
           IF(ITYPE(ICALC).EQ.161)THEN
            WRITE(*,*)'Reflecting atmospheric calculation'
           END IF
           NPATH1=ICALD(1,ICALC)
           NPATH2=ICALD(2,ICALC)
           I=1+NPATH2-NPATH1
           ICALATM=ICALD(4,ICALC)
           ICALCEL=ICALD(3,ICALC)
           WRITE(*,502)ICALATM,ICALCEL
502        FORMAT(' Atmosphere calculation : ',I3,
     1       ' Cell calculation : ',I3)
           WEIGHTF=BIT(2,ITYPE(ICALATM))
           EMISS=BIT(1,ITYPE(ICALATM))
           IF(WEIGHTF)THEN
            WRITE(*,*)'Weighting function'
           END IF
C          --------------------------------------------------------------
           IF(ITYPE(ICALCEL).EQ.128)THEN
             WRITE(*,*)'Single cell transmission'
             NSNGL=1 + ICALD(2,ICALC) - ICALD(1,ICALC)
             IF(WEIGHTF)THEN
               WRITE(*,501)NSNGL
501            FORMAT(' There are ',I3,' paths : ')
               DO J=NPATH1,NPATH2
                WRITE(*,*)J
               END DO
               CALL PROMPT(' Enter desired path : ')
               READ*,IPATH
               N=INT((1+NPATH2-NPATH1)/NSNGL)
               IPATH1=NPATH1 + (IPATH-1)*N
               IPATH2=IPATH1 + N-1
               GOTO 210
             ELSE
               WRITE(*,501)I
               CALL PROMPT(' Enter desired path : ')
               READ*,IPATH
               GOTO 131
             END IF
C          -------------------------------------------------------------
           ELSE IF(ITYPE(ICALCEL).EQ.129.OR.ITYPE(ICALCEL).EQ.131)THEN
             IF(ITYPE(ICALCEL).EQ.129)THEN
              WRITE(*,*)'PMR 2-pressure approximation'
             ELSE
              WRITE(*,*)'SCR'
             END IF
             CALL PROMPT(' Sideband(1) or Wideband(2) ? : ')
             READ*,J
             N=INT(0.5*(1+NPATH2-NPATH1))
             NPATH3=NPATH1 + N - 1
             IF(WEIGHTF)THEN
              IF(J.EQ.1)THEN
               IPATH1=NPATH1
               IPATH2=NPATH3
               CELLPATH = ICALD(1,ICALCEL)
              ELSE
               IPATH1=NPATH3+1
               IPATH2=NPATH2
               CELLPATH = ICALD(2,ICALCEL)
              END IF
              GOTO 210
             ELSE
              IF(J.EQ.1)THEN
 	       IPATH=NPATH1
              ELSE
               IPATH=NPATH2
              END IF
              GOTO 131
             END IF
C          --------------------------------------------------------------
           ELSE IF(ITYPE(ICALCEL).EQ.130)THEN
             WRITE(*,*)'PMR FILE'
             IPATH1=ICALD(1,ICALC)
             IPATH2=ICALD(2,ICALC)
             NPHAS=1+IPATH2-IPATH1
             WRITE(*,*)'Following combined atm/cell paths are
     1 calculated:'
             WRITE(*,*)'Phase angle    Path No.'
             DO I=1,NPHAS
              WRITE(*,*)(I-1)*360./NPHAS,IPATH1+I-1
             END DO
             WRITE(*,*)'Calculate individual path transmission (1) or '
             CALL PROMPT('combined SB and WB transmission (2) ? ')
             READ*,J
             IF(J.EQ.1)THEN
              CALL PROMPT('Enter path number : ')
              READ*,IPATH
              GOTO 131
             ELSE
              CALL PROMPT('Sideband(1) or Wideband(2) : ')
              READ*,IPMR
              GOTO 231
             END IF
           END IF
C        ------------------------------------------------------------------
        ELSE IF((ITYPE(ICALC).GE.128).AND.(ITYPE(ICALC).LE.131))THEN
          WRITE(*,*)'Cell calculation'
          IF(ITYPE(ICALC).EQ.128)THEN
          WRITE(*,*)'Single cell transmission'
          NSNGL=1 + ICALD(2,ICALC) - ICALD(1,ICALC)
          WRITE(*,501)NSNGL
          DO J=1,NSNGL
            I=ICALD(1,ICALC)+J-1
            WRITE(*,*)I
          END DO
          CALL PROMPT(' Enter desired path : ')
          READ*,IPATH
          EMISS=.FALSE.
          GOTO 131
C       ------------------------------------------------------------------
        ELSE IF(ITYPE(ICALC).EQ.129.OR.ITYPE(ICALC).EQ.131)THEN
           IF(ITYPE(ICALC).EQ.129)THEN
            WRITE(*,*)'PMR 2-pressure approximation'
           ELSE
            WRITE(*,*)'SCR'
           END IF
           CALL PROMPT(' Sideband(1) or Wideband(2) ? : ')
           READ*,J
           IF(J.EQ.1)THEN
            IPATH = ICALD(1,ICALC)
           ELSE
            IPATH = ICALD(2,ICALC)
           END IF
           EMISS=.FALSE.
           GOTO 131
C       ------------------------------------------------------------------
        ELSE IF(ITYPE(ICALC).EQ.130)THEN
           WRITE(*,*)'PMR FILE'
           FSTCEL=ICALD(1,ICALC)
           LSTCEL=ICALD(2,ICALC)
           NPHAS=1+LSTCEL-FSTCEL
           WRITE(*,*)'The following cell paths are calculated:'
           WRITE(*,*)'Phase angle    Path No.'
           DO I=1,NPHAS
            WRITE(*,*)(I-1)*360./NPHAS,FSTCEL+I-1
           END DO
           WRITE(*,*)'Calculate individual path transmission (1) or '
           CALL PROMPT('combined SB and WB transmission (2) ? ')
           READ*,J
           EMISS=.FALSE.
           IF(J.EQ.1)THEN
            CALL PROMPT('Enter path number : ')
            READ*,IPATH
            GOTO 131
           ELSE
            IPATH1=FSTCEL
            IPATH2=LSTCEL
            CALL PROMPT('Sideband(1) or Wideband(2) : ')
            READ*,IPMR
            GOTO 231
           ENDIF
           ENDIF
C        ------------------------------------------------------------------
        ELSE IF(ITYPE(ICALC).NE.200.AND.ITYPE(ICALC).NE.201)THEN
           WRITE(*,*)'Atmospheric calculation'
           IF(BIT(5,ITYPE(ICALC)))WRITE(*,*)'Reflecting paths'
           NPATH1=ICALD(1,ICALC)
           NPATH2=ICALD(2,ICALC)
           ANGLE=RCALD(1,ICALC)
           HT=RCALD(2,ICALC)
           WEIGHTF=BIT(2,ITYPE(ICALC))
           EMISS=BIT(1,ITYPE(ICALC))
           ABSORB=BIT(0,ITYPE(ICALC))
           IF(WEIGHTF)THEN
            WRITE(*,*)'Weighting function'
           END IF
           IF(EMISS)THEN
            WRITE(*,*)'Radiance calculation'
           END IF
           IF(ITYPE(ICALC).EQ.256)THEN
            WRITE(*,*)'Scattering calculation'
            EMISS=.TRUE.
            WEIGHTF=.FALSE.
           ENDIF
           IF(WEIGHTF)THEN
            CELLPATH=-1
            IPATH1=NPATH1
            IPATH2=NPATH2
            GOTO 210
           ELSE
            IF(NPATH2.EQ.NPATH1)THEN
             WRITE(*,*)'Plotting spectrum for PATH : ',NPATH1
             IPATH=NPATH1
            ELSE
             WRITE(*,*)'Paths run from ',NPATH1,' to ',NPATH2
             CALL PROMPT('Enter path : ')
             READ*,IPATH
            END IF
            GOTO 131
           END IF
        END IF

        IF(ITYPE(ICALC).EQ.200.OR.ITYPE(ICALC).EQ.201)THEN
         WRITE(*,*)'Reflecting layers calculation'
         IF(ITYPE(ICALC).EQ.201)THEN
          WRITE(*,*)'Combined atmosphere/cell paths'
         ELSE
          WRITE(*,*)'Pure atmosphere paths'
         END IF

         KE=ICALD(1,ICALC)
         KR=ICALD(2,ICALC)
         WRITE(*,366)KE
         WRITE(*,356)KR
366      FORMAT('Number of atmosphere calculations = ',I3)
356      FORMAT('Number of reflected atmosphere calcs = ',I3)
         WRITE(*,*)'Atmospheric calculations are:'
         DO K=1,KE
          WRITE(*,*)ICALD(2+K,ICALC)
         END DO
         WRITE(*,*)'Reflection Atmosphere calculations are:'
         DO K=1,KR
          WRITE(*,*)ICALD(2+K+KE,ICALC)
         END DO

         CALL PROMPT('How many reflecting calculations to include? : ')
         READ*,NREFL
         DO K=1,NREFL
          CALL PROMPT('Enter reflecting atm. calc. No. : ')
          READ*,IRATM
          CALL PROMPT('Enter albedo of layer : ')
          READ*,ALBEDO(K)
          NPATH1=ICALD(1,IRATM)
          NPATH2=ICALD(2,IRATM)
          SOLZEN(K) = RCALD(1,IRATM)
          IF(NPATH1.NE.NPATH2)THEN
           WRITE(*,*)'There is more than one path in this calculation'
           WRITE(*,397)NPATH1,NPATH2
397        FORMAT('NPATH1 = ',I3,' NPATH2 = ',I3)
           CALL PROMPT('Enter path number : ')
           READ*,IPATH
          ELSE
           IPATH=NPATH1
          END IF
          RPATH(K)=IPATH
         END DO

         CALL PROMPT('Enter distance from sun (AU) : ')
         READ*,SOLDIST
         DOMEGA=PI*(THETA0/SOLDIST)**2
         CALL PROMPT('Enter surface temperature of Sun (K) : ')
         READ*,TEMPSUN

         CALL PROMPT('How many other atm calcs to include? : ')
         READ*,NATMC
         DO K=1,NATMC
          CALL PROMPT('Enter atmosphere calc. No. : ')
          READ*,IATMC
          NPATH1=ICALD(1,IATMC)
          NPATH2=ICALD(2,IATMC)
          IF(NPATH1.NE.NPATH2)THEN
           WRITE(*,*)'There is more than one path in this calculation'
           WRITE(*,397)NPATH1,NPATH2
           CALL PROMPT('Enter path number : ')
           READ*,IPATH
          ELSE
           IPATH=NPATH1
          END IF
          APATH(K)=IPATH
         END DO


         GOTO 331

        END IF

       ELSE IF(Q.EQ.'Q')THEN
        STOP
       END IF
      GOTO 10
C---------------------------------------------------
C     section to plot a single spectrum
131   VMAX=VMIN+(NPOINT-1)*DELV
C      print*,VMAX,VMIN,EMISS
      COMPLE = .FALSE.
      IF(.NOT.EMISS)THEN
       CALL PROMPT('plot [1-output]?')
       READ(*,201)Q
       CALL UPCASE(Q)
       IF(Q.EQ.'Y')THEN
        COMPLE=.TRUE.
        ABSORB=.NOT.ABSORB
       END IF
      ELSE
       CALL PROMPT('plot brightness temperature? : ')
       READ(*,201)Q
       CALL UPCASE(Q)
       IF(Q.EQ.'Y')BRTEMP=.TRUE.
      END IF
      CALL PROMPT('Plot against wavenumber (1) or microns (2) : ')
      READ*,ISPACE
      IF(ISPACE.EQ.2)THEN
       MICRON=.TRUE.
      ELSE
       MICRON=.FALSE.
      END IF
      DO 110 I=1,NPOINT
      IREC=IREC0+(I-1)*NPATH+IPATH-1
C      print*,I,IPATH,NPATH,IREC0,IREC
      READ(2,REC=IREC)Y(I),ERROR(I)
C      print*,I,IREC,Y(I),ERROR(I)
      IF(COMPLE)Y(I)=1.-Y(I)
c      X(I)=VMIN+(I-1)*DELV
      X(I)=DVMIN+dble(I-1)*DDELV
      IF(BRTEMP)THEN
        IF(Y(I).GT.0)THEN
         Y(I)=BRIGHT(real(X(I)),Y(I))
        ELSE
         Y(I)=0.
        END IF
      ELSE
       IF(EMISS.AND.MICRON)Y(I)=Y(I)*real(X(I)*X(I))/1E4
      END IF
      IF(MICRON)X(I)=1E4/X(I)
110   CONTINUE
      CALL PROMPT('log plot? [Y/N]')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        LOG=.TRUE.
        YMAX=-50.0
        YMIN=50.0
       ELSE
        LOG=.FALSE.
        YMAX=Y(1)
        YMIN=YMAX
      END IF
      XMAX=real(X(1))
      XMIN=XMAX
      DO 304 I=1,NPOINT
      IF(LOG)Y(I)=PLOG10(Y(I))
      IF(Y(I).NE.-999.9)THEN
       YMAX=MAX(YMAX,Y(I))
       YMIN=MIN(YMIN,Y(I))
      ENDIF
      XMAX=MAX(XMAX,real(X(I)))
      XMIN=MIN(XMIN,real(X(I)))
304   CONTINUE
      IF(MICRON)THEN
       WRITE(*,184)XMIN,XMAX
      ELSE
       WRITE(*,204)XMIN,XMAX
      END IF
184   FORMAT(' wavelength ( x ) limits are: ',2F10.3)
204   FORMAT(' wavenumber ( x ) limits are: ',2F10.3)
      WRITE(*,305)YMIN,YMAX
305   FORMAT(' y limits are: ',2E12.5)
      IF(YMIN.EQ.YMAX)THEN
        Print*,'Ylimits are the same - no data'
        CALL PROMPT('input xmin,xmax,ymin,ymax')
        READ(*,*)XMIN,XMAX,YMIN,YMAX
      ELSE
        CALL PROMPT('force plotting axes? [y/n]')
        READ(*,201)Q
201     FORMAT(1A1)
        CALL UPCASE(Q)
        IF(Q.EQ.'Y')THEN
         CALL PROMPT('input xmin,xmax,ymin,ymax')
         READ(*,*)XMIN,XMAX,YMIN,YMAX
        ELSE
C        XMAX=VMAX
C        XMIN=VMIN
         YMAX=YMAX + (YMAX-YMIN)*0.05
        END IF
      END IF

C     *** start of debug lines ***
C      PRINT*,'ABSORB = ', ABSORB
C     *** end of debug lines ***

      IF(MICRON)THEN
       XLABEL='Wavelength [microns]'
      ELSE
       XLABEL='Wavenumbers [cm-1]'
      END IF
      IF(EMISS)THEN
       IF(BRTEMP)THEN
        YLABEL='Brightness Temperature [K]'
       IF(LOG)YLABEL='Log10(Brightness Temperature [K])'
       ELSE IF(MICRON)THEN
        YLABEL='Radiance [W cm-2 sr-1 micron-1]'
        IF(LOG)YLABEL='Log10(Radiance [W cm-2 sr-1 micron-1])'
       ELSE
        YLABEL='Radiance [W cm-2 sr-1 (cm-1)-1]'
        IF(LOG)YLABEL='Log10(Radiance [W cm-2 sr-1 (cm-1)-1])'
       END IF
      ELSE IF(ABSORB)THEN
       YLABEL='Absorption'
       IF(LOG)YLABEL='Log10(Absorption)'
      ELSE
       YLABEL='Transmission'
       IF(LOG)YLABEL='Log10(Transmission)'
      END IF
      IF(Q.EQ.'D')THEN
       YLABEL = 'Path spectrum'
       IF(LOG)YLABEL = 'Log10(Path spectrum)'
      END IF
      NPLOT=NPOINT
      GOTO 999
C------------------------------------------------------------------
C     this section to plot weighting functions and other complicated stuff

210   CONTINUE
C
C     first version just plots against log pressure
      VMAX=VMIN+(NPOINT-1)*DELV
      WRITE(*,211)NPOINT,VMIN,VMAX
211   FORMAT(' there are ',I8,' points covering the wavenumber range ',
     1F8.2,' - ',F8.2)
      CALL PROMPT('plot which wavenumber?')
      READ(*,*)VPLOT
      IPLOT=NINT((VPLOT-VMIN)/DELV)+1
      VPLOT=VMIN+(IPLOT-1)*DELV
      WRITE(*,213)IPLOT,VPLOT
213   FORMAT(' the nearest point is ',I8,' at ',F7.2,' cm-1')
      CALL PROMPT('Plot against height (1) or pressure(2) : ')
      READ(*,*)HTYPE
      IF(WEIGHTF)THEN
       N=1+IPATH2-IPATH1
       YMIN=1e30
       XMIN=1e30
       YMAX=-1e30
       XMAX=-1e30
       IF(CELLPATH.EQ.-1)THEN
        OLDOUT=1.
       ELSE
        IREC=IREC0 + (IPLOT-1)*NPATH + CELLPATH -1
        READ(2,REC=IREC)OLDOUT,ERROUT
       ENDIF
       I=0
       DO 233 J=IPATH1,IPATH2
        I=I+1
        IREC=IREC0+(IPLOT-1)*NPATH + J - 1
        READ(2,REC=IREC)OUTPUT,ERROUT
        X(I)=dble(OLDOUT-OUTPUT)
        IATP=LAYINC(NLAYIN(J),J)
        IF(HTYPE.EQ.1)THEN
          Y(I)=BASEH(IATP)
        ELSE
          Y(I)=BASEP(IATP)
        ENDIF
        OLDOUT=OUTPUT
        YMAX=MAX(YMAX,Y(I))
        YMIN=MIN(YMIN,Y(I))
233    CONTINUE
       IF(Y(1).EQ.Y(N))THEN
        NPLOT=INT(0.5*N)
        DO 224 I=1,NPLOT
         X(I)=X(I) + X(N+1-I)
         XMAX=MAX(XMAX,real(X(I)))
         XMIN=MIN(XMIN,real(X(I)))
224     CONTINUE
       ELSE
        NPLOT=N
        DO 332 I=1,NPLOT
         XMAX=MAX(XMAX,real(X(I)))
         XMIN=MIN(XMIN,real(X(I)))
332     CONTINUE
       END IF
       XLABEL='Contribution'
       IF(HTYPE.EQ.1)THEN
        YLABEL='Height [km] '
       ELSE
        YLABEL='Pressure [atm] '
       ENDIF
      ELSE
C      No other code ready yet
      END IF
      print*,'xmin,xmax,ymin,ymax'
      print*,xmin,xmax,ymin,ymax
      if(xmin.eq.xmax)then
       print*,'xmin = xmax. Enter new values'
       read*,xmin,xmax
      endif
      if(ymin.eq.ymax)then
       print*,'ymin = ymax. Enter new values'
       read*,ymin,ymax
      endif
      GOTO 999


231   CONTINUE
C     PMR cycles.

C     Cell transmission
      VMAX=VMIN+(NPOINT-1)*DELV
      DO 111 I=1,NPOINT

      WB=0.
      DO 113 J=1,NPHAS
       IPATH=IPATH1 + (J-1)
       IREC=IREC0+(I-1)*NPATH+IPATH-1
       READ(2,REC=IREC)Y1(J),ERROR(J)
       WB=WB+Y1(J)
113   CONTINUE
      WB=WB/NPHAS
      SB=0.
      DO 114 J=1,NPHAS
       Y1(J)=Y1(J)-WB
       IF(Y1(J).LT.0.)Y1(J)=-Y1(J)
       SB=SB+Y1(J)
114   CONTINUE
      SB=SB/NPHAS
      IF(IPMR.EQ.1)THEN
       Y(I)=SB
      ELSE
       Y(I)=WB
      END IF
      X(I)=DVMIN+(I-1)*DDELV
111   CONTINUE
      CALL PROMPT('log plot? [Y/N]')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        LOG=.TRUE.
        YMIN=50.0
        YMAX=-50.0
       ELSE
        LOG=.FALSE.
        YMIN=Y(1)
        YMAX=YMIN
      END IF
      XMIN=real(X(1))
      XMAX=real(X(1))
      DO 307 I=1,NPOINT
      IF(LOG)Y(I)=PLOG10(Y(I))
      IF(Y(I).NE.-999.9)THEN
       YMAX=MAX(YMAX,Y(I))
       YMIN=MIN(YMIN,Y(I))
      END IF
       XMAX=MAX(XMAX,real(X(I)))
       XMIN=MIN(XMIN,real(X(I)))
307   CONTINUE
      WRITE(*,234)VMIN,VMAX
234   FORMAT(' wavenumber ( x ) limits are: ',2F10.3)
      WRITE(*,325)YMIN,YMAX
325   FORMAT(' y limits are: ',2E12.5)
      CALL PROMPT('force plotting axes? [y/n]')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        CALL PROMPT('input xmin,xmax,ymin,ymax')
        READ(*,*)XMIN,XMAX,YMIN,YMAX
       ELSE
        XMAX=VMAX
        XMIN=VMIN
        YMAX=YMAX + (YMAX-YMIN)*0.05
        END IF
      XLABEL='Wavenumbers [cm-1]'
      YLABEL=' '
      NPLOT=NPOINT
      GOTO 999


331   CONTINUE
C     Reflecting layer models
      VMAX=VMIN+(NPOINT-1)*DELV
      CALL PROMPT('plot brightness temperature? : ')
      READ(*,201)Q
      CALL UPCASE(Q)
      BRTEMP=.FALSE.
      IF(Q.EQ.'Y')BRTEMP=.TRUE.
      CALL PROMPT('Plot against wavenumber (1) or microns (2) : ')
      READ*,ISPACE
      IF(ISPACE.EQ.2)THEN
       MICRON=.TRUE.
      ELSE
       MICRON=.FALSE.
      END IF
      DO 183 I=1,NPOINT
      X(I)=DVMIN+(I-1)*DDELV
      SFLUX=PLANCK(real(X(I)),TEMPSUN)*DOMEGA
      print*,x(i),'sflux : ',sflux
      Y(I)=0.

      DO K=1,IRATM
       IPATH=RPATH(K)
       IREC=IREC0+(I-1)*NPATH+IPATH-1
       READ(2,REC=IREC)YT,ERROUT
       Y(I)=Y(I)+SFLUX*YT*COS(SOLZEN(K)*PI/180.0)*ALBEDO(K)/PI
      END DO

      DO K=1,IATMC
       IPATH=APATH(K)
       IREC=IREC0+(I-1)*NPATH+IPATH-1
       READ(2,REC=IREC)YT,ERROUT
       Y(I)=Y(I)+YT
      END DO

      IF(BRTEMP)THEN
       IF(Y(I).GT.0)THEN
        Y(I)=BRIGHT(real(X(I)),Y(I))
       ELSE
        Y(I)=0.
       END IF
      END IF

C      print*,x(i),y(i)


      IF(MICRON)THEN
       Y(I)=Y(I)*real(X(I)*X(I))/1E4
       X(I)=1E4/X(I)
      END IF

C      print*,x(i),y(i)

183   CONTINUE
      CALL PROMPT('log plot? [Y/N]')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        LOG=.TRUE.
       ELSE
        LOG=.FALSE.
      END IF
      IF(LOG)THEN
       YMIN=PLOG10(Y(1))
      ELSE
       YMIN=Y(1)
      END IF
      YMAX=YMIN
      XMAX=real(X(1))
      XMIN=XMAX
      DO 374 I=1,NPOINT
      IF(LOG)Y(I)=PLOG10(Y(I))
      YMAX=MAX(YMAX,Y(I))
      YMIN=MIN(YMIN,Y(I))
      XMAX=MAX(XMAX,real(X(I)))
      XMIN=MIN(XMIN,real(X(I)))
374   CONTINUE
      IF(MICRON)THEN
       WRITE(*,184)XMIN,XMAX
      ELSE
       WRITE(*,204)XMIN,XMAX
      END IF
      WRITE(*,305)YMIN,YMAX
      CALL PROMPT('force plotting axes? [y/n]')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        CALL PROMPT('input xmin,xmax,ymin,ymax')
        READ(*,*)XMIN,XMAX,YMIN,YMAX
       ELSE
C        XMAX=VMAX
C        XMIN=VMIN
        YMAX=YMAX + (YMAX-YMIN)*0.05
        END IF
      IF(MICRON)THEN
       XLABEL='Wavelength [microns]'
      ELSE
       XLABEL='Wavenumbers [cm-1]'
      END IF
      IF(BRTEMP)THEN
        YLABEL='Brightness Temperature [K]'
        IF(LOG)YLABEL='LOG(Brightness Temperature [K])'
       ELSE IF(MICRON)THEN
        YLABEL='Radiance [W cm-2 sr-1 micron-1]'
        IF(LOG)YLABEL='LOG(Radiance [W cm-2 sr-1 micron-1])'
       ELSE
        YLABEL='Radiance [W cm-2 sr-1 (cm-1)-1]'
        IF(LOG)YLABEL='LOG(Radiance [W cm-2 sr-1 (cm-1)-1])'
      END IF
      NPLOT=NPOINT


      GOTO 999

999   CONTINUE


C-----------------------------------------------------------------------------
C     plotting section
C-----------------------------------------------------------------------------
	FIRST=.TRUE.

        OPEN(12,FILE='plspec.idl',STATUS='UNKNOWN')
C       print*,'nplot = ',nplot
        WRITE(12,*)NPLOT
        WRITE(12,*)'XLIM,YLIM : '
        WRITE(12,*)XMIN,XMAX,YMIN,YMAX
        WRITE(12,*)'XLABEL, YLABEL : '
        WRITE(12,322)XLABEL
        WRITE(12,322)YLABEL
322     FORMAT(A)
323     FORMAT(1X,F15.8,' ',E16.7)
        DO 998 I=1,NPLOT
                IF(real(X(I)).LT.XMIN)GOTO 998
                IF(real(X(I)).GT.XMAX)GOTO 998
C		WRITE(12,323)X(I),Y(I)
		WRITE(12,*)X(I),Y(I)
998     CONTINUE

        CLOSE(12)

	WRITE(*,*)' '
	WRITE(*,*)'***************** Pl_spec Complete *****************'

        STOP

        END

C------------------------------------------------------------------------------


      REAL FUNCTION PLOG10(Y)
      REAL Y

      IF(Y.GT.0)THEN
       PLOG10=LOG10(Y)
      ELSE
       PRINT*,'PLOG10. Trying to take log of ',Y
       PRINT*,'Setting to -999.9'
       PLOG10=-999.9
      ENDIF
      RETURN
      END
