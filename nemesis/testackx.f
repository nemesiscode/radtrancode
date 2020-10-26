      program testackx

      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE '../radtran/includes/constdef.f'
      character*100 ipfile,buffer
      integer imodel,iplanet,npro,nvmr,AMFORM,npro1,ncont
      real latitude,molwt,H(MAXPRO),P(MAXPRO),T(MAXPRO)
      REAL VMR(MAXPRO,MAXGAS),X,CONT(MAXCON,MAXPRO),flux,frain
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),I,J,JVMR
      REAL X1(MAXPRO),X2(MAXPRO),XMOL(MAXPRO),XVMR(MAXGAS)
      REAL CALCMOLWT,dist,albedo,Fsol,DENSCOND,RADCOND,MWCOND
      call reservegas

      call prompt('Enter name of input prf file : ')
      read(5,1)ipfile
1     format(a)


      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      MOLWT=-1.
C     First skip header
54    READ(1,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)AMFORM

      IF(AMFORM.EQ.0)THEN
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
      ELSE
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
      ENDIF

      IF(NPRO.GT.MAXPRO)THEN
          PRINT*,'Error. NPRO>MAXPRO ',NPRO,MAXPRO
          STOP
      ENDIF

      DO 20 I=1,NVMR
       READ(1,*)IDGAS(I),ISOGAS(I)
       PRINT*,I,IDGAS(I),ISOGAS(I)
20    CONTINUE

      CALL PROMPT('Enter gas number to apply scheme to (1-NVMR) : ')
      READ*,JVMR

C     Skip header
      READ(1,*)
      DO 30 I=1,NPRO
        READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
        IF(AMFORM.EQ.0)THEN
          XMOL(I)=MOLWT
        ELSE
          DO J=1,NVMR
           XVMR(J)=VMR(I,J)
          ENDDO
          XMOL(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
        ENDIF
30    CONTINUE

      CLOSE(UNIT=1)

      OPEN(UNIT=1,FILE='aerosol.prf',STATUS='OLD')

C     First skip header
55     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 55
       READ(BUFFER,*)NPRO1,NCONT
       IF (NPRO1.NE.NPRO)THEN
        print*,'NPRO <> NPRO1'
        stop
       ENDIF
       DO 31,I=1,NPRO
        READ(1,*)X,(CONT(J,I),J=1,NCONT)
31     CONTINUE
      CLOSE(1)


      call prompt('Enter imodel (0 = Lewis, 1=Ackerman) : ')
      READ(5,*)imodel

      call prompt('Enter distance of planet from sun (AU) and albedo: ')
      read(5,*)dist,albedo

      Fsol = albedo*1380/dist**2

      print*,'Incident solar flux (W/m2)  = ',Fsol
      print*,'Radiative equilibrium outgoing flux (W/m2)  = ',Fsol/4.0

      if(imodel.eq.1) then
       Call prompt('Enter flux (W m-2), frain : ')
       READ(5,*)flux,frain
      endif


      DENSCOND=1.0
      RADCOND=0.1
      MWCOND=12.0

      CALL ACKERMANMARLEYX1(IPLANET,LATITUDE,AMFORM,NPRO,NVMR,IDGAS,
     1 ISOGAS,P,T,H,VMR,XMOL,NCONT,CONT,FLUX,IMODEL,FRAIN,JVMR,
     2 DENSCOND,RADCOND,MWCOND,X1,X2)

      open(12,file='testackx.txt',status='unknown')
      write(12,*)NPRO
      write(12,*)flux,imodel,frain,jvmr
      do i=1,npro
       write(12,*)h(i),p(i),t(i),x1(i),x2(i)
      enddo
      close(12)
      end
