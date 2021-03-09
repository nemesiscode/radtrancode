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
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      call file(runname,ipfile,'prf')
      call rdmod(ipfile)

      ipfile='aerosol.prf'
      call rddmod(ipfile)

      if(idiag.gt.0)print*,'Check_profile called' 
      icheck=0
      do 10 i= 1,npro
       if(T(i).lt.0.0.or.isnan(T(i))) icheck=1
       sum=0.
       do 20 j=1,nvmr
        if (vmr(i,j).lt.0.0.or.isnan(vmr(i,j))) icheck=1
        sum=sum+vmr(i,j)
20     continue
       do 30 j=1,ncont
        if(dust(j,i).lt.0.0) then 
         dust(j,i)=0
         if(idiag.gt.0)then
          print*,'forcing dust to zero for icont,jlayer = ',j,i
         endif
        endif
        if(dust(j,i).lt.0.0.or.isnan(dust(j,i))) icheck=1
30     continue
10    continue

      if(idiag.gt.0)print*,'Checking : ',runname
      if(idiag.gt.0)print*,'icheck = ',icheck

      check_profile=icheck

      return

      end
