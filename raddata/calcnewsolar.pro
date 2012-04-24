openr,1,'solar.dat'
head=''
for i=1,3 do readf,1,head
data = fltarr(2,78)
readf,1,data
close,1
dist = 1.496e11
; Units of solar.dat are W m-2 nm-1. Want to convert to W um-1
data(1,*)=data(1,*)*4*!pi*dist^2  ; W nm-1
data(1,*)=data(1,*)*1e3		  ; W um-1

openw,1,'test.dat'
printf,1,data
close,1

end
