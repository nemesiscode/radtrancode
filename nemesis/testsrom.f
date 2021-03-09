      program testsrom
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      character*100 ipfile,head
      integer npro,i
      real patm(MAXPRO),temp(MAXPRO),xnew(MAXPRO)
      real PD,PT,RHC,RHM,VX,ch4tropvmr,ch4stratvmr
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      idiag=1

      print*,'Enter input test filename'
      read(5,1)ipfile
1     format(a)
   
      open(12,file=ipfile,status='old')
      read(12,*)npro
      do 10 i=1,npro
       read(12,*)patm(i),temp(i)
10    continue
      read(12,*)ch4tropvmr,ch4stratvmr
      read(12,*)PD,PT
      read(12,*)RHC,RHM
      read(12,*)VX
      close(12)

      call modifych4sromovsky(npro,patm,temp,ch4tropvmr,
     1  ch4stratvmr,PD,PT,RHC,RHM,VX,xnew)

      open(13,file='testsrom.txt',status='unknown')
      write(13,*)npro
      do 20 i=1,npro
       write(13,*)patm(i),temp(i),xnew(i)
20    continue
      close(13)

      if(idiag.gt.0)print*,'Run OK'
      end
