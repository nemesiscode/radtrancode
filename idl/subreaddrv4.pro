pro subreaddrv4,filename,nlay,npath,ngas,ncont,idgas,isogas,iprogas,$
baseh,delh,basep,baset,totam,press,temp,dopp,ppgas,amgas,cloud

openr,1,filename
head=' '
readf,1,head
readf,1,vmin,dv,npoint,fwhm
readf,1,x,y
readf,1,head
itmp=lonarr(4)
readf,1,itmp
imod=itmp(0)
ipara=itmp(1)
ncont=itmp(2)
flagc=itmp(3)
if(ncont gt 0) then readf,1,head
itmp=lonarr(3)
readf,1,itmp
nlay=itmp(0)
npath=itmp(1)
ngas=itmp(2)
idgas=intarr(ngas)
isogas=idgas
iprogas=idgas
for i=0,ngas-1 do begin
 readf,1,j
 idgas(i)=j

 readf,1,j,k
 isogas(i)=j
 iprogas(i)=k
endfor
for i=1,4 do readf,1,head
len = 9+ngas*2+ncont+ipara+flagc*(ncont+1)
data=fltarr(len,nlay)
readf,1,data
close,1

baseh = data(1,*)
delh = data(2,*)
basep = data(3,*)
baset = data(4,*)
totam = data(5,*)
press = data(6,*)
temp = data(7,*)
dopp = data(8,*)

ppgas = fltarr(ngas,nlay)
amgas = fltarr(ngas,nlay)

for i=0,ngas-1 do begin
 j = 9 + i*2
 for l=0,nlay-1 do begin
  amgas(i,l)=data(j,l)
  ppgas(i,l)=data(j+1,l)
 endfor
endfor

if (ncont gt 0) then begin
 cloud = fltarr(ncont,nlay)
 ic = j+2

 for i=0,ncont-1 do begin
  j = ic+i
  for l=0,nlay-1 do begin
   cloud(i,l)=data(j,l)
  endfor
 endfor

endif


return

end
