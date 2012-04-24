      program write_path
C     $Id:
C     ***********************************************************************
C
C     User friendly interactive set up of *.pat
C
C     Pat Irwin		25/10/94
C
C     ***********************************************************************
      character*100,infile,patfile,linedata,model,tmodel,filename
      character*100 dmodel,dspec
      character*15 scatter,therm,wf,cg,absorb,binbb,broad
      character*1,ans,ans1
      real layht,layang,angle,length,limit
      integer gid,pid,iid,laybot
      logical limb, pmr
      call prompt('Enter name of output .pat file : ')
      read(5,1)infile
      call file(infile,patfile,'pat')
    
      open(12,file=patfile,status='new')

      print*,'Instructions for using Path'
      print*,'============================='
      print*,' '
      print*,'Path processes a file full of ascii keywords (extension '
      print*,'.pat) and produces a "driver" file (extension .drv) for '
      print*,'use by genlbl. '
      print*,' ' 
      print*,'This program leads a new, or forgetful user, through the'
      print*,'input parameters to try and ease the use of Path.'
      print*,' '
      print*,'Press <RETURN> to start input procedure ...'
      print*,' '
      read(5,2)ans
      print*,'Enter desired model : '
      print*,'(1) LBL'
      print*,'(2) BAND'
      print*,'(3) CORRK'
      call prompt('Enter choice : ')
      read*,imod
      write(12,*)' '
      write(12,202)'interval'
      call prompt('Enter wavenumber min, max and output spacing and
     1 resolution : ')
      read*,vmin,vmax,delv,fwhm
      write(12,*) vmin,vmax,delv,fwhm
      print*,' '
      if(imod.eq.1)then
C     LBL Model
       if(fwhm.eq.0.)then
        iconv=0
       else
        iconv=1
       end if
       call prompt('Enter wing limit and overall limit : ')
       read*,wing,limit
      else if(imod.eq.2)then
C      Band Model
       print*,'Select model'
       print*,'------------'
       print*,'1) Overlapping lines - (line data)'
       print*,'2) General random band model - (line data)'
       print*,'3) Malkmus-Lorentz model - (band data)'
       print*,'4) Goody-Lorentz model - (band data)'
       print*,'5) Godson-Lorentz model - (band data)'
       print*,'6) Goody-Voigt model - (band data)'
       read*,i
       iconv=9+i
       if(i.lt.3)then
         print*,' '
         print*,'Equivalent Width Calculation'
         print*,'Enter code for combining equivalent widths'
         print*,'0 - Equivalent width used is Lorentz only'
         print*,'1 - W=sqrt(WL^2 + WD^2 - (WL*WD/Sm)^2)'
         print*,'2 - Interpolate from Voigt EW Look-up Table'
         call prompt('Enter combination code : ')
         read*,com_mod
         print*,' '
         if(com_mod.ne.2)then
          print*,'Enter code for calculating Lorentz widths'
          print*,'(meaningless if previous choice=2)'
          print*,'0 - Equiv. Lorentz = WEAK'
          print*,'1 - Equiv. Lorentz = STRONG'
          print*,'2 - Equiv. Lorentz = COMBINED (Rodgers and 
     1Williams 1974)'
          call prompt('Enter Lorentz code : ')
          read*,lor_mod
         else
          lor_mod = 0.
         end if
         if(com_mod.eq.1)then
          print*,'Enter code for calculating Doppler widths'
          print*,'0 - Equiv. Doppler = WEAK'
          print*,'1 - Equiv. Doppler = STRONG'
          print*,'2 - Equiv.Doppler = COMBINED (Rodgers and 
     1Williams 1974)'
          call prompt('Enter Doppler code : ')
          read*,dop_mod
         else
          dop_mod=0
         end if

         wing = dop_mod + lor_mod*3 + com_mod*9
         limit = i-1
       else
         wing=0.
         limit=0.
       end if
      else
C      Corr_K Model
       wing=0.
       limit=0.
      end if
      write(12,*) iconv,wing,limit

      write(12,202)' '
202   format(a)

      print*,'Enter specdata filename (linedata or band data'
      call prompt('depending on model choice) : ')
      read(5,1)linedata
1     format(a60)
      write(12,202)'spec data'
      write(12,202)linedata
      write(12,202)' '
      print*,' '

      call prompt('Enter an atmosphere ? : ')
      read(5,2)ans
      if(ans.eq.'Y'.or.ans.eq.'y')then

      print*,'Enter name of atmospheric temperature and vmr profiles'
      call prompt('file (prf extention may be omitted) : ')
      read(5,1)model
      print*,' '

      tmodel='model  '
      tmodel(7:)=model(1:)
      write(12,202) tmodel
      write(12,202)' '

      call prompt('Enter dust profile? (Y/N) : ')
      read(5,2)ans
      if(ans.eq.'Y'.or.ans.eq.'y')then
       call prompt('Enter name of dust profile : ')
       read(5,1)model
       dmodel='dust model  '
       dmodel(12:)=model(1:)
       write(12,202)dmodel
       call prompt('Enter name of dust spectral file : ')
       read(5,1)model
       dmodel='dust spectra  '
       dmodel(14:)=model(1:)
       write(12,202)dmodel
       write(12,202)' '
      end if

      print*,'Layer Specification'
      print*,'========================='
      print*,'path splits the atmosphere into a number of layers. This'
      print*,'is separate from calculations of the atmospheric paths so'
      print*,'that you can have several paths through a set of layers'
      print*,'(eg for multiple tangent heights).'
      print*,' '

      write(12,202)'layer'

      call prompt('Enter no. of layers to split the atmosphere into :')
      read*,nlay
      write(12,203)nlay
203   format('nlay ',i3)
      print*,' '

      call prompt('Enter the height of the base of the lowest layer :')
      read*,layht
      write(12,204)layht
204   format('layht ',f8.2)
      print*,' '

      print*,'Enter the angle from the zenith for splitting'
      print*,'into layers. For example you would use zero for nadir'
      print*,'calculations but 90 for limb. the stored layers are'
      print*,'scaled by the cos of this angle so that they are always'
      print*,'similar to the nadir case. You might use layang=0 for'
      print*,'multiple limb paths or for simulating someone elses'
      call prompt('code. Enter angle : ')
      read*,layang
      write(12,205)layang
205   format('layang ',f5.1)
      print*,' '

      print*,'Enter an integer specifying the layering type.'
      print*,'0=split by equal changes in pressure over layers'
      print*,'1=split by equal changes in log pressure over layers'
      print*,'2=split by equal changes in height over layers'
      print*,'3=split by equal changes in path length at layang'
      call prompt('Enter layer code : ')
      read*,laytyp
      write(12,206)laytyp
206   format('laytyp ',i3)
      print*,' '

      print*,'Enter an integer specifying how to integrate'
      print*,'vertically over the layer. '
      print*,'0=use values at centre of layer,' 
      print*,'1=use curtis godson paths for a gas with constant mixing'
      print*,'  ratio and integrate using Simpsons rule'
      print*,'1 is usually better but 0 useful for intercomparisons.'
      call prompt('Enter integration code : ')
      read*,layint
      write(12,207)layint
207   format('layint ',i3)
      print*,' '
      print*,'Input of layers complete'
      print*,' '

      print*,'Do you want to perform multiple atmospheric '
      call prompt('calculations ? : ')
      read(5,2)ans
      if(ans.eq.'Y'.or.ans.eq.'y')then

      print*,'Atmospheric Path Specification'
      print*,'=============================='

      print*,' '
      call prompt('Enter number of paths through the atmosphere : ')
      read*,natm
      print*,' '
      do 10 i=1,natm

       print*,'Atmospheric Path : ',i
     
       write(12,202)' '
       write(12,202)'atm'
       call prompt('Limb(1) or Nadir(2) : ')
       read*,j
       if(j.eq.1)then
        limb=.true.
        call prompt('Enter number of bottom layer : ')
        read*,laybot
        write(12,208)laybot
208     format('limb ',i3)
       else
        limb=.false.
        call prompt('Enter number of bottom layer : ')
        read*,laybot
        call prompt('Enter angle from zenith : ')
        read*,angle
        write(12,209)angle,laybot
209     format('nadir ',f5.1,i3)
       end if

       call prompt('Weighting Function Calculation (Y/N)? : ')
       read(5,2)ans
2      format(a1)
       if(ans.eq.'y'.or.ans.eq.'Y')then
        wf='wf'
       else
        wf='nowf'
       end if

       if(imod.ne.2)then
       call prompt('Use Curtis Godson Approximation (Y/N)? : ')
       read(5,2)ans
       if(ans.eq.'y'.or.ans.eq.'Y')then
        cg='cg'
       else
        cg='nocg'
       end if
       else
        print*,'Curtis-Godson paths required by Band models'
        cg='cg'
       end if

       therm='notherm'

       if(wf.eq.'wf')then
        scatter='noscatter'
        therm='notherm'
        absorb='noabsorb'
        broad='nobroad'
        binbb='nobinbb'
       else
        call prompt('Thermal Emission Calculation (Y/N)? : ')
        read(5,2)ans

        if(ans.eq.'y'.or.ans.eq.'Y')then
         therm='therm'
         print*,'Calculate emission (i.e. Planck function) outside'
         print*,'Radtran to simulate broad band calculations? '
         call prompt('(Y/N) : ')
         read(5,2)ans
         if(ans.eq.'y'.or.ans.eq.'Y')then
          broad='broad'
         else
          broad='nobroad'
         end if
         print*,'Calculate planck function at bin centres not each'
         call prompt('wavenumber - much faster (Y/N)? : ')
         read(5,2)ans
         if(ans.eq.'y'.or.ans.eq.'Y')then
          binbb='binbb'
         else
          binbb='nobinbb'
         end if
         absorb='noabsorb'
         scatter='noscatter'
        else 
         call prompt('Scattering calculation ? ')
         read(5,2)ans
         if(ans.eq.'y'.or.ans.eq.'Y')then
          scatter='scatter'
          therm='notherm'
          broad='nobroad'
          binbb='nobinbb'
          absorb='noabsorb'
          cg='nocg'
         else
          scatter='noscatter'
          therm='notherm'
          broad='nobroad'
          binbb='nobinbb'
          print*,'Calculate absorption rather than transmission? Needed'
          print*,'for small absorptions because of the use of adaptive'
          call prompt('integration. (Y/N)? : ')
          read(5,2)ans
          if(ans.eq.'y'.or.ans.eq.'Y')then
           absorb='absorb'
          else
           absorb='noabsorb'
          end if
         end if
         end if
       end if
 
       write(12,202)therm
       write(12,202)scatter
       write(12,202)wf
       write(12,202)cg
       write(12,202)absorb
       write(12,202)binbb
       write(12,202)broad

10     continue

      end if

      if(cg.eq.'nocg')then
      call prompt('Add a gas cell (Y/N)? : ')
      read(5,2)ans
      if(ans.eq.'y'.or.ans.eq.'Y')then
        write(12,202)'cell'
        call prompt('Enter number of gases in cell : ')
        read*,ngas
        write(12,210)ngas
210     format('gases ',i3)

        print*,'Enter identifier, isotope and volume mixing ratio of' 
        print*,'each gas. The identifier is NOT always the same as '
        print*,'Hitran. Isotope zero includes all isotopes at '
        print*,'terrestrial ratios. Genlbl corrects line strengths for'
        print*,'explicitly specified isotopes as if they were 100%'
        do 20 j=1,ngas
         print*,'gas ',j
         call prompt('Enter ID, ISO and VMR : ')
         read*,id,iso,vmr
         write(12,*)id,iso,vmr
20      continue
        call prompt('Enter Dopler shift for cell (cm-1) : ')
        read*,dop
        write(12,211)dop
211     format('dop ',e12.5)
        call prompt('Enter cell length (cm) : ')
        read*,length
        write(12,212)length
212     format('length ',e12.5)

        call prompt('Single path cells(1) or PMC(2) : ')
        read*,j
        if(j.eq.1)then
         write(12,202)'sngl'
         call prompt('Enter number of single paths : ')
         read*,npath
         write(12,*)npath
         do 24 k=1,npath
          print*,'path ',k
          call prompt('Enter temperature(K) and pressure(atm) : ')
          read*,temp,press
          write(12,*)temp,press
24       continue
        else
         pmr=.true.
         call prompt('Use two pressure approx.(1) or file(2) : ')
         read*,k
         if(k.eq.1)then
          write(12,202)'pmr twop'
          call prompt('Enter low temperature(K) and pressure(atm) : ')
          read*,temp,press
          write(12,*)temp,press
          call prompt('Enter high temperature(K) and pressure(atm) : ')
          read*,temp,press
          write(12,*)temp,press
         else
          write(12,202)'pmr file'
          call prompt('Enter PMC cycle filename : ')
          read(5,1)filename
          write(12,202)filename
         end if
        end if
       end if

      end if
      end if

      print*,' '
      call prompt('Enter a single atmospheric layer ? : ')
      read(5,2)ans
      if(ans.eq.'y'.or.ans.eq.'Y')then

       print*,' '
       print*,'Defining single atmospheric layers'
       print*,'=================================='

       print*,'Genlbl treats atmospheric and cell layers differently so'
       print*,'this is not the same as defining single cell paths. '
       print*,'The input details may either be read in from an external'
       print*,'file or included in the .pat file. '
       write(12,202)' '
       call prompt('Specify a filename (Y/N) : ')
       read(5,2)ans
       if (ans.eq.'y'.or.ans.eq.'Y')then
        call prompt('Enter filename : ')
        read(5,1)filename
        tmodel='sngatm  '
        tmodel(8:)=filename
        write(12,202)tmodel
       else
        write(12,202)'sngatm'
        print*,'Enter parameters explicitly'
        print*,' '
        call prompt('Calculate emission (Y/N)? : ')
        read(5,2)ans
        if(ans.eq.'y'.or.ans.eq.'Y')then
         absorb='therm'
        else
         print*,'Calculate absorption rather than transmission? Needed'
         print*,'for small absorptions because of the use of adaptive'
         call prompt('integration. (Y/N)? : ')
         read(5,2)ans1
         if(ans1.eq.'y'.or.ans1.eq.'Y')then
          absorb='absorb'
         else
          absorb='noabsorb'
         end if
        end if
        write(12,202)absorb

        print*,'Enter pressure(atm), temperature(K) and number of'
        call prompt('gases : ')
        read*,press,temp,ngas
        write(12,*)ngas,press,temp

        print*,'The first line in the file should contain the keyword'
        print*,'absorb or noabsorb. The following line contains the '
        print*,'pressure, temperature and number of gases. Each gas is'
        print*,'For each gas enter the identifier, isotope,'
        print*,'volume mixing ratio and amount. If the amount'
        print*,'is < 1.e10 it is assumed to be a path length in km.'

        do 99 i=1,ngas
         print*,'gas no. ',i
         call prompt('Enter id,iso,vmr and amount : ')
         read*,id,iso,vmr,amount
         write(12,*)id,iso,vmr,amount
99      continue

        call prompt('Enter number of dust types : ')
        read*,ncont
        write(12,*)ncont
        if(ncont.gt.0)then
         call prompt('Enter dust x-section filename : ')
         read(5,202)dmodel
         write(12,202)dmodel
         do i=1,ncont
          print*,'Dust type ',i
          call prompt('Enter dust amount (number/cm2) : ')
          read*,amount
          write(12,*)amount
         end do
        end if

       end if

      end if


      if(imod.eq.1)then
       print*,' '
       call prompt('Enter integ. error limit (Y/N)? (Default 1%) : ')
       read(5,2)ans
       if(ans.eq.'y'.or.ans.eq.'Y')then
        call prompt('Enter error : ')
        read*,error
        write(12,220)error
220     format('error ',e12.5)
       end if     

       print*,'Define how to process line data ? '
       print*,'================================='

       print*,'If chosen, user enters the gas id, isotope id and an'
       print*,'integer to define how to do the processing. The default'
       print*,'is zero and uses voigt lines etc. These parameters are'
       print*,' hacked into genlbl as needed and force funny line'
       print*,'wings etc.'

       call prompt('Enter process parameters (Y/N)? : ')
       read(5,2)ans
       if(ans.eq.'y'.or.ans.eq.'Y')then
        call prompt('Enter gas id, isotope id and process integer : ')
        read*,gid,iid,pid
        write(12,221)gid,iid,pid
221     format('process ',i3,i3,i3)
       end if     

      end if

      call prompt('Remove unused layers (Y/N)? : ')
      read(5,1)ans
      write(12,202)' '
      if(ans.eq.'y'.or.ans.eq.'Y')then
        write(12,202)'clrlay'
      else
        write(12,202)'noclrlay'
      end if

      if(imod.ne.2)then

       call prompt('Combine all cell/atm paths (Y/N)? : ')
       read(5,1)ans
       write(12,202)' '
       if(ans.eq.'y'.or.ans.eq.'Y')then
        write(12,202)'combine'
       end if

      else
        write(12,202)'nocombine'
      end if


      end

