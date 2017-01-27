openr,1,'partfextra.dat'
openw,2,'partfextra.new'

head=''
for i=0,1 do begin
 readf,1,head
 printf,2,head
endfor
ngas=1
readf,1,ngas
printf,2,ngas
idiso=intarr(2,ngas)
readf,1,idiso
printf,2,idiso

ans=''
for igas=0,ngas-1 do begin
  for i=0,1 do begin
   readf,1,head
   printf,2,head
  endfor
  nl=1
  readf,1,nl
  printf,2,nl
  tmp=dblarr(2,nl)
  readf,1,tmp
  printf,2,tmp
  print,head
  plot_io,tmp(0,*),tmp(1,*)
  read,ans
endfor

close,/all

end
