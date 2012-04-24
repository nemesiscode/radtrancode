      subroutine check_limits(nlayer,npath,ngas,ncont,nsec,
     1 ng,npk,ntk,nlayin,layinc)
C     **************************************************************
C     Subroutine to check that array sizes are not exceeded.
C
C     Pat Irwin	30/7/01	Original
C     Pat Irwin 29/2/12 Updated for Radtrans2.0
C
C     **************************************************************
      implicit none
      INCLUDE '../includes/arrdef.f'
      integer nlayer,npath,ngas,ncont,nsec, ng,npk,ntk
      integer nlayin(maxpat),layinc(maxlay,maxpat),ipath,i

      if (nlayer.gt.maxlay) then
                write (*,*) ' CHECK_LIMITS: Too many layers'
                write (*,*) ' Nlayer = ',nlayer,' maxlay = ',maxlay
                stop
      endif

      if (npath.gt.maxpat) then
                write (*,*) ' CHECK_LIMITS: Too many paths'
                write (*,*) ' Npath = ',npath,' maxpat = ',maxpat
                stop
      endif

      if (ngas.gt.maxgas) then
                write (*,*) ' CHECK_LIMITS: Too many gases'
                write (*,*) ' Ngas = ',ngas,' Maxgas = ',maxgas
                stop
      endif
                
      if (ncont.gt.maxcon) then
                write (*,*) ' CHECK_LIMITS: Too many dust continua'
                write (*,*) ' Ncont = ',ncont,' maxcon = ',maxcon
                stop
      endif

      if (nsec.gt.maxsec) then
                write (*,*) ' CHECK_LIMITS: Too many dust continua pts'
                write (*,*) ' Nsec = ',nsec,' Maxsec = ',maxsec
                stop
      endif

      if (ng.gt.maxg) then
                write (*,*) ' CHECK_LIMITS: Too many g ordinates'
                write (*,*) ' Ng = ',ng,' Maxg = ',maxg
                stop
      endif

      if ((npk.gt.maxk).or.(ntk.gt.maxk)) then
                write (*,*) ' CHECK_LIMITS: Too many P/T points in K',
     1                  ' tables'
                write (*,*) ' Npk = ',npk,' Maxk = ',maxk
                write (*,*) ' Ntk = ',ntk,' Maxk = ',maxk
                stop
      endif
                
      do Ipath = 1, npath
                if (nlayin(Ipath).gt.maxinc) then
                  write (*,*) ' CHECK_LIMITS: Nlayin exceeds maxinc'
                  write (*,*) ' Ipath = ', Ipath,' Nlayin = ',
     1                          nlayin(Ipath), ' maxinc = ', maxinc
                  stop
                endif
                do I = 1, nlayin(Ipath)
                        if (layinc(I,Ipath).gt.nlayer) then
                          write (*,*) ' CHECK_LIMITS: Layinc exceeds',
     1                                ' nlayer'
                          write (*,*) ' Ipath = ', Ipath,
     1                                  ' Layer = ', I, ' Layinc = ',
     2                                  layinc(I,Ipath), ' Nlayer = ',
     3                                  nlayer
                          stop
                        endif
                enddo
      enddo
 
      return

      end
