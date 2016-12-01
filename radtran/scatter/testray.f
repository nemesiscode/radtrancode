      program testray
      real v,T,P,x1,x2,x3,rayleigha,rayleighv,rayleighj

      call prompt('Enter wavenumber, pressure and temperature : ')
      read*,v,T,P
  
      x1=rayleighj(v,P,T)
      x2=rayleighv(v,P,T)
      x3=rayleigha(v,P,T)

      print*,'RayleighJ = ',x1
      print*,'RayleighV = ',x2
      print*,'RayleighA = ',x3
    
      end
