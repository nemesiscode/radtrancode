
      subroutine fpread(vv,temp,kpara,absh2,abshe,dah2dT,dahedT,
     1 deltaf)
C     **********************************************************************
C     Reads in the collision-induced absorption coefficients of molecular
C     hydrogen and helium for a range of different para-H2 fractions.
C
C     Units of absorption : cm-1
C
C     input parameters:
C     
C	vv            	real	is the wavenumber wanted (may not
C                              	be an exact value on the ciatable,
C                              	hence interpolate.
C
C     	temp         	real	the temperature, also may require
C                        	interpolation
C
C
C     Output Parameters 
C
C  	kpara		integer	Number if para-H2 fractions listed.
C
C	absh2(NPARA)	real	H2-H2 CIA absorption for different para-H2
C				fractions
C	abshe(NPARA)	real	H2-He CIA absorption for different para-H2
C				fractions
C	dah2dT(NPARA)	real	Rate of change of H2-H2 absorption with
C				temperature
C	dahedT(NPARA)	real	Rate of change of H2-He absorption with
C				temperature
C	deltaf		real	Step in para-H2 table
C
C     P.Irwin  (26/7/01)  Original
C     Fletcher (07/10/10) Removed explicit reference to 12 para fractions,
C				Changed to NPARA (also needs to be modified in 
C				nparacon.f and ciatable/fparatable.f
c     Fletcher (31/01/11) Added parameter DF to subroutine fpread.f

C     **********************************************************************

      implicit none
      include '../includes/ciacom.f'

      real deltaf
C     input/output variables
      real vv,temp,absh2(NPARA),abshe(NPARA),dah2dT(NPARA),dahedT(NPARA)
      REAL KH2H21(12,NUMT,NUMWN),KH2HE1(12,NUMT,NUMWN)
C     Other variables
      real fp1(12)

      real abs11, abs12, abs21, abs22
      real bbs11, bbs12, bbs21, bbs22
      real u,t,temp1,dudt
      integer cv,ct,pindex,i,j,k
      integer kpara
      logical fexist
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet
C     ----------------------------------------------------------------------

C     define the fixed parameters

      kpara=ipara2

      if(iread2.eq.1)then
        call datarchive(aname2)
      
        inquire(file=aname2,exist=fexist)
        if(fexist)then
         if(idiag.gt.0)print*,'fpread. Reading: ',aname2
        else
         if(idiag.gt.0)print*,'fpread. File does not exist : ',aname2
        endif
    
        open(12,file=aname2,form='unformatted',status='old')      
         if(idiag.gt.0)print*,'FPREAD. Opening look-up file'
         if(idiag.gt.0)print*,'File =',aname2
         if(idiag.gt.0)print*,'Tabulated temperatures : '
         read(12)tempk2
         do 21 i=1, NUMT
            if(idiag.gt.0)print*,i,tempk2(i)
21       continue
         if(idiag.gt.0)print*,'Tabulated para-H2 fractions : '
C        Older para-h2 tables have 12 ortho/para ratios instead of 12
         if(ipara2.eq.12)then
          read(12)fp1
          do i=1,ipara2
           fp(i)=fp1(i)
          enddo
         else
          read(12)fp
         endif

         do 22 i=1, ipara2
            if(idiag.gt.0)print*,i,fp(i)
22       continue
       
         if(ipara2.eq.12)then
          read(12)kh2h21
          read(12)kh2he1
          do 301 i=1,ipara2
           do 302 j=1,NUMT
            do 303 k=1,NUMWN
             kh2h2(i,j,k)=kh2h21(i,j,k)
             kh2he(i,j,k)=kh2he1(i,j,k)
303         continue
302        continue
301       continue
         else
          read(12)kh2h2
          read(12)kh2he
         endif
         close(12)

C     ----------------------------------------------------------------------
C     initialisation section
C     ----------------------------------------------------------------------

C        set up the wavenumber array
         do 23 i=1,NUMWN
            vvk2(i) = (i-1)*DNU2
23       continue

C        initialise the absorption array to be zero in case of failure
         do 24 i=1, IPARA2
            absh2(i) = 0.0
            abshe(i) = 0.0
   	    dah2dt(i)=0.0
	    dahedt(i)=0.0
 24      continue

         if(idiag.gt.0)print*,'fpread : New Data read in'

         iread2=0

      end if

C     ----------------------------------------------------------------------
C     range checking
C     ----------------------------------------------------------------------

      if(vv.gt.vvk2(NUMWN).or.vv.lt.vvk2(1).and.ivk.ne.1)then
         if(idiag.gt.0)print*,'fpread : wavenumber out of range'
         ivk=1
         do i=1,npara
          absh2(i)=0.0
          abshe(i)=0.0
 	  dah2dt(i)=0.0
	  dahedt(i)=0.0
         enddo
         return
      end if

      if(temp.gt.tempk2(NUMT).and.it1k.ne.1)then
         if(idiag.gt.0)print*,'fpread: warning temp > MaxT '
         if(idiag.gt.0)print*,'Requested Temperature = ',temp
         if(idiag.gt.0)print*,'Maxumum temperature = ',tempk2(NUMT)
         temp1 = tempk2(NUMT)
         it1k=1
      else
         temp1 = temp
      end if

      if(temp.lt.tempk2(1).and.it2k.ne.1)then
         if(idiag.gt.0)print*,'fpread: warning temp <  MinT ',temp
         if(idiag.gt.0)print*,'Setting to ',tempk2(1)
         temp1 = tempk2(1)
         it2k=1
      else
         temp1 = temp
      end if

C     ----------------------------------------------------------------------
C     main section
C     ----------------------------------------------------------------------

      deltaf=fp(2)-fp(1)

    
C     go away and find the plave in the table for wavenumber:
      i = NUMWN
      cv = pindex(vvk2,i,vv)

C     ... and for temperature:
      i = NUMT
      ct = pindex(tempk2,i,temp1)

C     The returned values are the indices of the array immediately
C     below or equal to the required temp or wavenumber. These indices
C     are the same, of course, for both H2-H2 and H2-He gas pairs.

      t = (vv - vvk2(cv))/(vvk2(cv+1)-vvk2(cv))
      u = (temp1 - tempk2(ct))/(tempk2(ct+1)-tempk2(ct))
      dudt = 1.0/(tempk2(ct+1)-tempk2(ct))


C     Loop over all tabulated para-H2 fractions
      do 26 i=1,IPARA2
         abs11 = kh2h2(i,ct,cv)
         abs12 = kh2h2(i,ct,cv+1)
         abs22 = kh2h2(i,ct+1,cv+1)
         abs21 = kh2h2(i,ct+1,cv)
         bbs11 = kh2he(i,ct,cv)
         bbs12 = kh2he(i,ct,cv+1)
         bbs22 = kh2he(i,ct+1,cv+1)
         bbs21 = kh2he(i,ct+1,cv)
 
         absh2(i) = (1.0-t)*(1.0-u)*abs11 + t*(1.0-u)*abs12 + 
     &        t*u*abs22 + (1.0-t)*u*abs21
         abshe(i) = (1.0-t)*(1.0-u)*bbs11 + t*(1.0-u)*bbs12 + 
     &        t*u*bbs22 + (1.0-t)*u*bbs21
	 dah2dT(i) = -(1.0-t)*abs11 -t*abs12 + t*abs22 + (1.0-t)*abs21
	 dahedT(i) = -(1.0-t)*bbs11 -t*bbs12 + t*bbs22 + (1.0-t)*bbs21

	 dah2dt(i)=dah2dt(i)*dudt
	 dahedt(i)=dahedt(i)*dudt

26    continue


      return

      end

