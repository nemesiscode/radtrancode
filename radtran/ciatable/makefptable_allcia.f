      program makefptable_allcia
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
      include '../includes/ciacom.f'

C k is the master-table of absorption co-efficients
      character*100 opfile

      dimension f(601),alf(601), temps(NUMT)
      double precision temps,temp
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      idiag=1

      print*,'Program MAKEFPTABLE_ALLCIA - calculates pressure-induced '
      print*,'        absorption co-efficients for ortho/para H2-H2 and'
      print*,'        H2-He, plus H2-N2, N2-CH4, N2-N2, CH4-CH4, H2-CH4'
      print*,' '


      print*,'Enter Tmin and Tstep for calculation'
      read*,TMIN,TSTEP

      print*,'Enter number if fpara fractions (12 or 24)'
      read*,NUMPARA

      FMIN = 0.25
      print*,'Enter DF for para-table'
      read*,dfin
      
      DF = dfin*FLOAT(NPARA)/FLOAT(NUMPARA)
C      DF = 0.01*FLOAT(NPARA)/FLOAT(NUMPARA) 

      print*,'Enter dnu (usually 1.0 for fparatables)'
      read*,DNU

      print*, 'Enter width of slit for N2-H2, N2-N2 and H2-CH4: '
      print*, '[Need 4.4 cm-1 or larger]: '
      read*,slit


      print*,'Enter name of output file : '
      read(5,1)opfile
1     format(a)

      print*, 'Temperature range is: ', TMIN, ' to ',
     &         TMIN+(TSTEP*NUMT)




CCCCCCCCCCCCC
CCCCCCCCCCCCC Calculation of Borysow Ortho-Para CIA Fractions
CCCCCCCCCCCCC

C initialise the table of absorption-co-efficients by setting
C all entries = 0

       do 6 i1=1,NUMPARA
        do 7 i2=1,NUMT
         do 8 i3=1,NUMWN
          kh2h2(i1,i2,i3)= 0.0
          kh2he(i1,i2,i3)= 0.0
8        continue
7       continue
6      continue

       do 200 itemp = 1,NUMT
	  TEMP = TMIN + (itemp - 1)*TSTEP
          TEMPK2(ITEMP) = SNGL(TEMP)

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

 201    print*,'Ortho/Para calculation complete'


CCCCCCCCCCCCC
CCCCCCCCCCCCC Calculation of Borysow CIA for other gases
CCCCCCCCCCCCC

C five gas pairs at present:
C
C 1............H2-N2
C 2............N2-CH4
C 3............N2-N2
C 4............CH4-CH4
C 5............H2-CH4



C initialise the table of absorption-co-efficients by setting
C all entries = 0
       do 16 i1=1,NUMPAIRS-4
        do 17 i2=1,NUMT
         do 18 i3=1,NUMWN
          kcia5(i1,i2,i3)= 0.0
18       continue
17      continue
16     continue


       do 300 itemp = 1,NUMT
	 temps(itemp) = TMIN + (itemp - 1)*TSTEP
	  tempk1(itemp)=SNGL(temps(itemp))

	  print*,'temp = ',temps(itemp)

C Fletcher 04/2018:  In below N2-H2, N2-CH4, etc, had to multiply by DNU to get 
C Correct coefficients from the subroutines.


C =========================== N2 - H2 ===========================
C Used n2h2_sub (Teanby) rather than n2h2_s which had some memory issues
 355	  igas = 1
	  temp=temps(itemp)
	  do istep = 1,3
	   fnumin=0.0 + (istep-1)*500.0*dnu
	   fnumax = fnumin + 499.0*dnu
	   nf = 0
	   call n2h2_sub(temp,fnumin,fnumax,dnu,nf,f,alf,slit)
	   do i=1,nf
	      j=1 + int(f(i)/dnu)
	      kcia5(igas,itemp,j)=kcia5(igas,itemp,j)+sngl(alf(i))
	   end do
	  enddo
	  print*,'n2h2 OK'


C =========================== N2 - CH4 ===========================
 356	  igas = 2
	  temp=temps(itemp)
	  do istep = 1,2
	   fnumin=0.0 + (istep-1)*500.0*dnu
	   fnumax = fnumin+499.0*dnu
	   nf = 0
	   call n2ch4_s(temp,fnumin,fnumax,dnu,nf,f,alf)
	   do i=1,nf
	      j=1 + int(f(i)/dnu)
	      kcia5(igas,itemp,j)=kcia5(igas,itemp,j)+sngl(alf(i))
	   end do
	  enddo
	  print*,'n2ch4 OK'

C =========================== N2 - N2 ===========================
 357	  igas = 3
	  temp=temps(itemp)
	  do istep = 1,2
	   fnumin=0 + (istep-1)*500.0*dnu
	   fnumax = fnumin + 499.0*dnu
	   nf = 0
	   call n2n2_s(temp,fnumin,fnumax,dnu,nf,f,alf,slit)
	   do i=1,nf
	      j=1 + int(f(i)/dnu)
	      kcia5(igas,itemp,j)=kcia5(igas,itemp,j)+sngl(alf(i))
	   end do
	  enddo
	  print*,'n2n2 OK'

C =========================== CH4 - CH4 ===========================
 358	  igas = 4
	  temp=temps(itemp)
	  do istep = 1,3
	   fnumin=0 + (istep-1)*500.0*dnu
	   fnumax = fnumin + 499.0*dnu
	   nf = 0
	   call ch4ch4_s(temp,fnumin,fnumax,dnu,nf,f,alf)
	   do i=1,nf
	      j=1 + int(f(i)/dnu)
	      kcia5(igas,itemp,j)=kcia5(igas,itemp,j)+sngl(alf(i))
	   end do
	  enddo
	  print*,'ch4ch4 OK'

C =========================== CH4 - H2 ===========================
 359	  igas = 5

C NB ********* there was a problem with h2ch4_s.f
C the CH4 octopole spectrum is going wrong at T > 150 K (approx)
C The effect gets worse with temp (try plotting it - strange spikes
C are appearing). Anyway, leaving this code in for the time being
C because the error is seemingly small.
C Problematic program replaced with h2ch4_sub. Output no agrees
c with Borysows original and ouput doesn't contain spikes. NT 23/6/04

          temp=temps(itemp)
          do istep = 1,3
           fnumin=0 + (istep-1)*500.0*dnu
           fnumax = fnumin + 499.0*dnu
           nf = 0
           call h2ch4_sub(temp,fnumin,fnumax,dnu,nf,f,alf,slit)
           do i=1,nf
              j=1 + int(f(i)/dnu)
              
              kcia5(igas,itemp,j)=kcia5(igas,itemp,j)+sngl(alf(i))
           end do
          enddo
          print*,'h2ch4 OK'


 300    continue

 301    print*,'Extra Gas Calculation complete'



      open(12,file=opfile,form='unformatted',status='unknown')

      write(12) tempk2
      write(12) fp
      write(12) kh2h2
      write(12) kh2he
      write(12) kcia5
      close(12)

      end


