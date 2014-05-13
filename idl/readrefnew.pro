pro readrefnew,afile,iform,nplanet,xlat,npro,ngas,molwt,idiso,height,$
press,temp,vmr

openr,1,afile
lhead:
head=''
readf,1,head
a=strmid(head,0,1)
if(a eq '#') then goto,lhead

iform=0
reads,head,iform
nsamp=1
readf,1,nsamp
tmp = fltarr(5)
if(iform eq 1) then tmp=fltarr(4)
readf,1,tmp
nplanet=long(tmp(0))
xlat=tmp(1)
npro = long(tmp(2))
ngas = long(tmp(3))
molwt=-999.9
if(iform eq 0) then molwt = tmp(4)
idiso = lonarr(2,ngas)
readf,1,idiso
head=''
readf,1,head

n1 = 3+ngas

tmp = fltarr(n1,npro)
readf,1,tmp
height=fltarr(npro)
press=height
temp=height
vmr=fltarr(npro,ngas)
height(*)=tmp(0,*)
press(*)=tmp(1,*)
temp(*)=tmp(2,*)
for i=0,ngas-1 do vmr(*,i)=tmp(i+3,*)

close,1

return

end



