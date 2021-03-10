      program makeciatable
C ====================== Rationale ===============================
C
C Program to create a look-up table of pressure-induced absorption
C co-efficients at a range of temperatures and wavenumbers for
C H2 and He gas-pairs. Program runs a series of sub-programs
C to calculate the opacities for each gas. Programs by Borysow et al.
C
C created P. Irwin 1/96. AOPP, University of Oxford, England.
C Pat Irwin	2/3/12	Checked for Radtrans2.0
C
C ================ Modification History ==========================
C
C modified to include N2-H2, N2-N2, CH4-CH4, CH4-H2 and N2-CH4
C pair opacities by  C.Nixon, 1/97. AOPP, Oxford.
C
C (note: UNIX you must compile all the .f files together into one
C executable. Subprograms:
C
C h2h2_v1s.f  h2he_v0s.f h2h2_v0s.f  h2h2_v2s.f  h2he_v1s.f
C n2n2_s.f n2h2_s.f n2ch4_s.f ch4ch4_s.f h2ch4_s.f
C
C ================================================================

C variable declarations
      implicit double precision (a-h,o-z)
      include '../includes/ciacom.f'

C k is the master-table of absorption co-efficients
      
C f is the frequency array, alf is the absorption co-efficient array
C and temps is the temperature array
      dimension f(601),alf(601), temps(NUMT)
      character*100 opfile
      double precision temps,temp
C nine gas pairs at present:
C
C 1............H2-H2 (ortho:para = 1:1 `equilibrium')
C 2............H2-He                        "
C 3............H2-H2 (ortho:para = 3:1 `normal')
C 4............H2-He                        "
C 5............H2-N2
C 6............N2-CH4
C 7............N2-N2
C 8............CH4-CH4
C 9............H2-CH4
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      idiag=1

      print*,'Program CIATABLE - calculates pressure-induced'
      print*,'                   absorption co-efficients.'
      print*,' '


      print*,'Enter Tmin and Tstep for calculation'

      print*,'Setups that are currently used'
      print*,'Borysowsub.tab uses Tmin=50, Tstep=15.'
      print*,'NewBorysow.tab uses Tmin=50, Tstep=15.'
      print*,'exociatable.tab uses Tmin=400, Tstep=100'
      print*,'exociatable_fletcher.tab uses Tmin=400, Tstep=100'

      read*,TMIN,TSTEP
    
      print*,'Enter dnu (usually 1.0 for longwave ISO and 10 otherwise)'
      read*,DNU


      print*,'Enter name of output file : '
      read(5,1)opfile
1     format(a)

      print*, 'Temperature range is: ', TMIN, ' to ',
     &         TMIN+(TSTEP*NUMT)

      print*, 'Enter width of slit for N2-H2, N2-N2 and H2-CH4: '
      read (5,*) slit

C initialise the table of absorption-co-efficients by setting
C all entries = 0
       do 6 i1=1,NUMPAIRS
        do 7 i2=1,NUMT
         do 8 i3=1,NUMWN
          kcia(i1,i2,i3)= 0.0
8        continue
7       continue
6      continue

       do 200 itemp = 1,NUMT
	  temps(itemp) = TMIN + (itemp - 1)*TSTEP
          tempk1(itemp)=SNGL(temps(itemp))

          print*,'temp = ',temps(itemp)

C this loop for equilibrium and normal ortho:para H2
          do 250, normal=0,1

             if (normal .eq. 0) then
                igas = 1
             else
                igas = 3
             endif
             print*, 'normal =', normal

             temp=temps(itemp)
             fnumin=0.
             fnumax=4000.
             nf = 0
             call h2h2_v0s(normal,temp,fnumin,fnumax,dnu,nf,f,alf)
             print*,'h2h2_v0 OK'

             do i=1,nf
                j=1 + int(f(i)/DNU)
                kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
             end do

             temp=temps(itemp)
             fnumin=2000.
             fnumax=7990
             call h2h2_v1s(normal,temp,fnumin,fnumax,dnu,nf,f,alf)
             print*,'h2h2_v1 OK'
             do i=1,nf
                j=1 + int(f(i)/dnu)
                kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
             end do

             temp=temps(itemp)
             fnumin=7000.
             fnumax=13000.
             call h2h2_v2s(normal,temp,fnumin,fnumax,dnu,nf,f,alf)
             print*,'h2h2_v2 OK'
             print*,'nf = ',nf
             do i=1,nf
                j=1 + int(f(i)/dnu)
                kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
                print*,f(i),j,alf(i)
             end do
             stop

             if (normal .eq. 0) then
                igas = 2
             else
                igas = 4
             endif

             temp=temps(itemp)
             fnumin=0.
             fnumax=4000.
             call h2he_v0s(normal,temp,fnumin,fnumax,dnu,nf,f,alf)
             print*,'h2he_v0 OK'
             do i=1,nf
                j=1 + int(f(i)/dnu)
                kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
             end do

             temp=temps(itemp)
             fnumin=2000.
             fnumax=7990
             call h2he_v1s(normal,temp,fnumin,fnumax,dnu,nf,f,alf)
             print*,'h2he_v1 OK'
             do i=1,nf
                j=1 + int(f(i)/dnu)
                kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
             end do

 250      continue

C =========================== N2 - H2 ===========================
 255      igas = 5
          temp=temps(itemp)
          fnumin=0
          fnumax = 1500
          nf = 0
C          print*,'Calling n2h2_s. slit = ',slit
C          call n2h2_s(temp,fnumin,fnumax,dnu,nf,f,alf,slit)
C          print*,'n2h2 OK'
C          do i=1,nf
C             j=1 + int(f(i)/dnu)
C             kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
C          end do

C =========================== N2 - CH4 ===========================
 256      igas = 6
          temp=temps(itemp)
          fnumin=0
          fnumax = 1000
          nf = 0
          call n2ch4_s(temp,fnumin,fnumax,dnu,nf,f,alf)
          print*,'n2ch4 OK'
          do i=1,nf
             j=1 + int(f(i)/dnu)
             kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
          end do

C =========================== N2 - N2 ===========================
 257      igas = 7
          temp=temps(itemp)
          fnumin=0
          fnumax = 1000
          nf = 0
          call n2n2_s(temp,fnumin,fnumax,dnu,nf,f,alf,slit)
          print*,'n2n2 OK'
          do i=1,nf
             j=1 + int(f(i)/dnu)
             kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
          end do
          temp=temps(itemp)

C =========================== CH4 - CH4 ===========================
 258      igas = 8
          temp=temps(itemp)
          fnumin=0
          fnumax = 1500
          nf = 0
          call ch4ch4_s(temp,fnumin,fnumax,dnu,nf,f,alf)
          print*,'ch4ch4 OK'
          do i=1,nf
             j=1 + int(f(i)/dnu)
             kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
          end do

C =========================== CH4 - H2 ===========================
 259      igas = 9

C NB ********* there is a problem with this code:
C the CH4 octopole spectrum is going wrong at T > 150 K (approx)
C The effect gets worse with temp (try plotting it - strange spikes
C are appearing). Anyway, leaving this code in for the time being
C because the error is seemingly small.

          temp=temps(itemp)
          fnumin=0
          fnumax = 1500
          nf = 0
          call h2ch4_s(temp,fnumin,fnumax,dnu,nf,f,alf,slit)
          print*,'h2ch4 OK'
          do i=1,nf
             j=1 + int(f(i)/dnu)
             kcia(igas,itemp,j)=kcia(igas,itemp,j)+sngl(alf(i))
          end do

 200    continue

 201    print*,'Calculation complete'

      open(12,file=opfile,form='unformatted',status='unknown')

      write(12) temps
      write(12) kcia
      close(12)

      end


