      subroutine ciaread(vv,temp,abs,dabsdT)
C     **********************************************************************
C     Reads in a collision-induced absorption coefficient table
C
C     Units of absorption : cm-1
C
C     Pat Irwin		Original	9/6/95
C			Revision	31/1/12
C
C
C     input parameters:
C
C     vv                 real       is the wavenumber wanted (may not
C                                   be an exact value on the ciatable,
C                                   hence interpolate.
C
C     temp               real       the temperature, also may require
C
C     Output variables
C       abs(9)          real    Absorption coefficients for 9-pairs
C                               (listed in nciacon.f)
C       dabsdT(9)       real    Rate of change of abs with temperature
C                                   interpolation
C     
C     **********************************************************************

      implicit none
      
C     define the fixed parameters
      include '../includes/ciacom.f'

C     input/output variables
      real vv,temp,abs(NUMPAIRS),dabsdT(NUMPAIRS)

C     The temp. array, the wavenumber array, the abs. co-eff array
      double precision temps(NUMT)
C     Note that array sizes must be hard-wired and must agree with the 
C     sizes printed out by ciatable.f.
C     The actual temperature values are variable (read-in).

      real abs11, abs12, abs21, abs22
      real u,t,temp1,dudt
      integer cv,ct,pindex,i,j,ii,jj,kk
      logical fexist

C     ----------------------------------------------------------------------

C     Elements of absorption table are:
C       1............H2-H2 (ortho:para = 1:1 `equilibrium')
C       2............H2-He                        "
C       3............H2-H2 (ortho:para = 3:1 `normal')
C       4............H2-He                        "
C       5............H2-N2
C       6............N2-CH4
C       7............N2-N2
C       8............CH4-CH4
C       9............H2-CH4


      if(iread1.eq.1)then

         call datarchive(aname1)

         inquire(file=aname1,exist=fexist)

         if(fexist)then
           print*,'ciaread. Reading File = ',aname1
         else
           print*,'Ciaread. File does not exist',aname1
           stop
         endif
           

         open(12,file=aname1,form='unformatted',status='old')

           read(12)temps
C          convert from double precision to real
           print*,'CIAREAD. Temperatures'
           do 21 i=1, NUMT
             tempk1(i) = sngl(temps(i))
             print*,tempk1(i)
21         continue

           read(12)kcia

         close(12)

C        ----------------------------------------------------------------------
C        initialisation section
C        ----------------------------------------------------------------------

C        this sets up the wavenumber array
         do 22 i=1,NUMWN
            vvk1(i) = (i-1)*dnu1
22       continue

C        initialise the absorption array to be zero in case of failure
         do 24 i=1, NUMPAIRS
            abs(i) = 0.0
            dabsdT(i)=0.0
24       continue
         print*,'ciaread : New Data read in'

C        Set read flag to zero to prevent table being read again
         iread1=0

      end if

C     ----------------------------------------------------------------------
C     range checking
C     ----------------------------------------------------------------------

      if((vv.gt.vvk1(NUMWN).or.vv.lt.vvk1(1)).and.ivc.ne.1)then
         print*,'ciaread : wavenumber out of range'
         ivc=1
         return
      end if

      if(temp.gt.tempk1(NUMT).and.it1c.ne.1)then
         print*,'ciaread: warning temp > MaxT ',temp
         print*,'Setting to ',tempk1(NUMT)
         temp1 = tempk1(NUMT)
         it1c=1
      else
         temp1 = temp
      end if

      if(temp.lt.tempk1(1).and.it2c.ne.1)then
         print*,'ciaread: warning temp <  MinT ',temp
         print*,'Setting to ',tempk1(1)
         temp1 = tempk1(1)
         it2c=1
      else
         temp1 = temp
      end if

C     ----------------------------------------------------------------------
C     main section
C     ----------------------------------------------------------------------

C     go away and find the place in the table ...
C     ... for wavenumber ...
      i = NUMWN
      cv = pindex(vvk1,i,vv)
C     ... and for temperature:
      i = NUMT
      ct = pindex(tempk1,i,temp1)

C     The returned values are the indices of the array immediately
C     below or equal to the required temp or wavenumber. These indices
C     are the same, of course, for all the NUMPAIRS gas pairs.

C     Loop over each gas pair and linearly interpolate the actual
C     absorption co-efficient

      do 26 i=1, NUMPAIRS
         abs11 = kcia(i,ct,cv)
         abs12 = kcia(i,ct,cv+1)
         abs22 = kcia(i,ct+1,cv+1)
         abs21 = kcia(i,ct+1,cv)

         t = (vv - vvk1(cv))/(vvk1(cv+1)-vvk1(cv))
         u = (temp1 - tempk1(ct))/(tempk1(ct+1)-tempk1(ct))
         dudt = 1.0/(tempk1(ct+1)-tempk1(ct))

         abs(i) = (1.0-t)*(1.0-u)*abs11 + t*(1.0-u)*abs12 +
     &        t*u*abs22 + (1.0-t)*u*abs21

         dabsdT(i) = (-(1.0-t)*abs11 -t*abs12 +
     &        t*abs22+(1.0-t)*abs21)*dudt

 26   continue

      return

      end

