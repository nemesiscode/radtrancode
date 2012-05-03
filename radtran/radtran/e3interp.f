      subroutine e3interp(xtrans,e3out,de3out)
C     ********************************************************
C     Subroutine for first reading table of E3 vs transmission and
C     then interpolating for the required transmission. Routine also 
C     returns the rate of change of E3 with transmission.
C
C     Input variable:
C	xtrans	double precision	Transmission to calculate E3
C
C     Output variables:
C	e3out   double precision	Interpolated E3
C	de3out	double precision	Interpolatd rate of change of 
C					E3 with transmission
C
C     Pat Irwin		27/9/10
C
C     ********************************************************
      implicit none
      integer mtab,n,i,j
      parameter (mtab=1001)
      double precision xtrans,e3out,de3out,tabtran(mtab),tabe3(mtab)
      double precision tau,tran,e3,egrad(mtab),dtrans
      double precision f        

      character*100 ipfile,buffer
      common /e3table/tabtran,tabe3,egrad

      dtrans=1.0/float(mtab-1)

      if(tabtran(mtab).ne.1.0) then
       print*,'Reading e3 table'
       ipfile='e3tab.dat'
       call datarchive(ipfile)
       open(12,file=ipfile,status='old')
       read(12,*)n
       read(12,1)buffer
1      format(a)
       do i=1,n
        read(12,*)tau,tran,e3
        tabtran(i)=tran
        tabe3(i)=e3
       enddo
       close(12)

       egrad(1)=(tabe3(2)-tabe3(1))/dtrans
       egrad(mtab)=(tabe3(mtab)-tabe3(mtab-1))/dtrans
       do i=2,mtab-1
        egrad(i)=0.5*(tabe3(i+1)-tabe3(i-1))/dtrans
       enddo

      endif
    

      j=1+int(xtrans/dtrans)
      if(j.ge.mtab)j=mtab-1


      f=(xtrans-tabtran(j))/dtrans
      if(f.gt.1.0)f=1.0
C      print*,xtrans,j,tabtran(j),f

      e3out = (1.0-f)*tabe3(j)+f*tabe3(j+1)
      de3out = (1.0-f)*egrad(j)+f*egrad(j+1)

      return

      end

