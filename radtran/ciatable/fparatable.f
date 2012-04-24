      program fparatable
C ====================== Rationale ===============================
C
C Program to create a look-up table of pressure-induced absorption
C co-efficients at a range of temperatures and wavenumbers for
C H2 and He gas-pairs for a range of para fractions:
C
C created P. Irwin 7/01. AOPP, University of Oxford, England.
C Pat Irwin	2/3/12	Checked for Radtrans2.0
C
C ================ Modification History ==========================
C
C (note: UNIX you must compile all the .f files together into one
C executable. Subprograms:
C
C nh2h2.f  nh2he.f
C
C ================================================================

C variable declarations
      implicit double precision (a-h,o-z)

C k is the master-table of absorption co-efficients
      real kh2h2(12,25,1501)
      real kh2he(12,25,1501)
      character*100 opfile
      real temps(25),fp(12)

C f is the frequency array, alf is the absorption co-efficient array
C and temps is the temperature array
      NUMPARA=12

C set limits of calculation grid and step size
      TMIN=50.0          ! you can change this
      TSTEP=15.0         ! and this, but not
      NUMTEMPS=25        ! this.

      NUMWN=1501         ! don't change this
      DNU=1.0            ! or this, because then other programs
C                                    which read the table will get confused

      FMIN = 0.25
      DF = 0.02 


      print*,'Program FPARATABLE - calculates pressure-induced'
      print*,'                   absorption co-efficients.'
      print*,' '
      print*,'Enter name of output file : '
      read(5,1)opfile
1     format(a)

      print*, 'Temperature range is: ', TMIN, ' to ',
     &         TMIN+(TSTEP*NUMTEMPS)


C initialise the table of absorption-co-efficients by setting
C all entries = 0

       do 6 i1=1,NUMPARA
        do 7 i2=1,NUMTEMPS
         do 8 i3=1,NUMWN
          kh2h2(i1,i2,i3)= 0.0
          kh2he(i1,i2,i3)= 0.0
8        continue
7       continue
6      continue

       do 200 itemp = 1,NUMTEMPS
	  temps(itemp) = TMIN + (itemp - 1)*TSTEP
          print*,'temp = ',temps(itemp)
          TEMP = TEMPS(ITEMP)

C this loop for equilibrium and normal ortho:para H2
          do 250 IPARA=1,NUMPARA
            FPARA = FMIN + (IPARA-1)*DF
            IF(IPARA.EQ.NUMPARA)FPARA=-1.0D0	! Equilibrium

            fp(IPARA)=SNGL(FPARA)
            print*,'fp = ',fp(IPARA)
            do 260 iwave = 1,NUMWN
             FREQ = DNU*(IWAVE-1)

             CALL NH2H2(FPARA,TEMP,FREQ,ALFA1)
             CALL NH2HE(FPARA,TEMP,FREQ,ALFA2)

             kh2h2(ipara,itemp,iwave)=sngl(ALFA1)
             kh2he(ipara,itemp,iwave)=sngl(ALFA2)
 260        continue

 250      continue

 200    continue

 201    print*,'Calculation complete'

      open(12,file=opfile,form='unformatted',status='unknown')

      write(12) temps
      write(12) fp
      write(12) kh2h2
      write(12) kh2he

      close(12)

      end


