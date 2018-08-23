      integer function check_profile(runname)
C     ****************************************************************
C     Code to read in latest prf file and see if it has non-negative 
C     temperatures and vmrs.
C
C     ****************************************************************

      implicit none
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE '../radtran/includes/pathcom.f'
C ../includes/pathcom.f holds the variables used by the software when
C calculating atmospheric paths (e.g. NPATH and IMOD).
      INCLUDE '../radtran/includes/laycom.f'
      integer i,j,icheck
      real sum
      character*100 runname,ipfile

      call file(runname,ipfile,'prf')
      call rdmod(ipfile)

      ipfile='aerosol.prf'
      call rddmod(ipfile)

      print*,'Check_profile called' 
      icheck=0
      do 10 i= 1,npro
       if(T(i).lt.0.0.or.isnan(T(i))) icheck=1
       print*,T(i),(vmr(i,j),j=1,nvmr)
       sum=0.
       do 20 j=1,nvmr
        if (vmr(i,j).lt.0.0.or.isnan(vmr(i,j))) icheck=1
        sum=sum+vmr(i,j)
20     continue
       print*,'Summed vmrs = ',sum
       print*,(dust(j,i),j=1,ncont)
       do 30 j=1,ncont
        if(dust(j,i).lt.0.0) then 
         dust(j,i)=0
         print*,'forcing dust to zero for icont,jlayer = ',j,i
        endif
        if(dust(j,i).lt.0.0.or.isnan(dust(j,i))) icheck=1
30     continue
10    continue

      print*,'Checking : ',runname
      print*,'icheck = ',icheck

      check_profile=icheck

      return

      end
