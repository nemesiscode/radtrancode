      program testack
      character*100 ipfile,opfile,aefile,qcfile
      real flux,Teff,frain
      integer imodel

      call prompt('Enter name of input prf file : ')
      read(5,1)ipfile
1     format(a)

      opfile='testack.prf'
      aefile='testaer.prf'
      qcfile='testqc.prf'

      Teff = 124
      flux = 5.67e-8*Teff**4.

      call prompt('Enter imodel, f_rain : ')
      READ(5,*)imodel,frain

      call ackermanmarley(ipfile,opfile,aefile,qcfile,flux,
     1  imodel,frain)

      end
