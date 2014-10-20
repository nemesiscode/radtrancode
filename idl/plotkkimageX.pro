; ********************************************************************
; IDL procedure to plot the dr/dx spectra contained in a standard kk.out
; file dumped by Nemesis retrieval code.
;
; Pat Irwin	28/10/03
;
; ********************************************************************
mx = 401
my = 1024
filename=' '
print,'Enter name of kk.out file : '
;filename='kk.out'
read,filename
openr,1,filename,/f77_unformatted
y = fltarr(my)
yn = y
kk = fltarr(my,mx)

readu,1,y,yn
readu,1,kk

close,1

print,'Enter runname'
runname=''
read,runname

filename = strcompress(string(runname,'.spx'),/remove_all)
openr,1,filename
xtemp = fltarr(4)
fwhm=0.0
readf,1,xtemp
fwhm = xtemp(0)
latitude = xtemp(1)
longitude = xtemp(2)
ngeom = long(xtemp(3))
nconv=1
nav = lonarr(ngeom)
for igeom = 0,ngeom-1 do begin
  readf,1,nconv
  if(igeom eq 0)then begin
   data=fltarr(ngeom,3,nconv)
   data1 = fltarr(3,nconv)
   vconv = fltarr(nconv)
  endif
  i1=1
  readf,1,i1
  nav(igeom)=i1
  tang = fltarr(6,i1)
  readf,1,tang
  readf,1,data1
  for j=0,2 do data(igeom,j,*)=data1(j,*)
endfor
ny = ngeom*nconv

close,1

print,'ny = ',ny

vconv = fltarr(nconv)
vconv(*)=data1(0,*)

print,'Range of vconv : ',data1(0,0),data1(0,nconv-1)

print,'Enter desired output xrange : '
a=1.0
b=1.0
read,a,b

sum=0
!p.position = 0
!p.multi=[0,1,ngeom]

for igeom=0,ngeom-1 do begin
 c = 0.0
 ioff = nconv*igeom
 for i=0,nconv-1 do begin
  if(vconv(i) ge a and vconv(i) le b ) then begin
   c = max([c,data(igeom,1,i),yn(i+ioff)])
  endif
 endfor

 plot,vconv,data(igeom,1,*),xrange=[a,b],yrange=[0,c],$
 xtitle='Wavenumber (cm-1)',ytitle='Radiance W cm-2 sr-1 cm'
 oplot,data(igeom,0,*),yn(ioff:(ioff+nconv-1)),linestyle=1

 
 for i=0,nconv-1 do sum = sum+(data(igeom,1,i) - yn(i+ioff))^2

endfor

sum = sqrt(sum/ny)
print,'rms error = ',sum


ans=' '
read,ans

print,'Print kk to window (0), postscipt(1), eps(2)'
ips = 0
read,ips

head=''
filename = strcompress(string(runname,'.prf'),/remove_all)
openr,1,filename
tmp = fltarr(5)
iform=0
readf,1,iform
if(iform eq 1) then tmp=fltarr(4)
readf,1,tmp
npro = long(tmp(2))
ngas = long(tmp(3))
for i=1,ngas+1 do readf,1,head
tmp = fltarr(3+ngas,npro)
readf,1,tmp
close,1

tref  = fltarr(npro)
pref  = fltarr(npro)

tref = tmp(2,*)
pref = tmp(1,*)

filename = strcompress(string(runname,'.apr'),/remove_all)
openr,1,filename
head=' '
readf,1,head
nvar=1
readf,1,nvar
varident = lonarr(nvar,3)
vartmp = lonarr(3)
aprdat=fltarr(nvar,3,npro)
for ivar=0,nvar-1 do begin
 readf,1,vartmp
 print,vartmp
 varident(ivar,*)=vartmp(*)
 if(vartmp(0) lt 500 and vartmp(2) eq 0) then begin
  ipfile=''
  readf,1,ipfile
  openr,2,ipfile
    npro = 1
    clen=1.5
    readf,2,npro,clen
    aprdat1 = fltarr(3,npro)
    readf,2,aprdat1
    for j=0,2 do aprdat(ivar,j,*)=aprdat1(j,*)
  close,2
 endif
 if(vartmp(2) eq 1)then for i=0,2 do readf,1,head
 if(vartmp(2) eq 2)then readf,1,head
 if(vartmp(2) eq 3)then readf,1,head
 if(vartmp(2) eq 9)then for i=1,3 do readf,1,head
 if(vartmp(2) eq 10)then for i=1,3 do readf,1,head
 if(vartmp(0) eq 555)then readf,1,head
 if(vartmp(0) eq 999)then readf,1,head

endfor
close,1


press = fltarr(npro)
press(*)=1.013*pref(0,*)

ioff=0
print,'Enter pressure range (max,min) : '
pmin = 0.0
pmax = 0.0
read,pmax,pmin
lp1 = alog(pmin)
lp2 = alog(pmax)

nl = 500
lpref = lp1+(lp2-lp1)*findgen(nl)/float(nl-1)
lpress = alog(reverse(press))

dv = vconv(1)-vconv(0)
vmin = vconv(0)
vmax = vconv(nconv-1)
nout = 1+long((vmax-vmin)/dv)
bb = fltarr(nout,nl)

print,'Number of variables are : ',nvar
for i=0,nvar-1 do print,varident(i,*)

print,'Enter nplotx, nploty'
nplotx = 2
nploty = 2
read,nplotx,nploty

print,'Enter colour table (33 is a good one!) : '
ict=0
read,ict

iplotx = 1
iploty = 1

if(ips ne 0) then begin
 set_plot,'ps'
 if(ips eq 1) then begin
  device,filename='plotkkimage.ps',encapsulated=0
 endif else begin
  device,filename='plotkkimage.eps',/encapsulated
 endelse

 xwid = 75
 ywid = 45
; DEVICE,/COLOR,BITS=8,/PORT, xsize=xwid, ysize=ywid,xoffset=2,$
;                yoffset=1
 DEVICE,/COLOR,BITS=8,/PORT, xsize=xwid, ysize=ywid
        !p.position=0  
        !p.multi=0

endif else set_plot,'x'

for ivar = 0,nvar-1 do begin
 print,'Variable : ',ivar,varident(ivar,*)
 print,'ioff = ',ioff
 vartyp = varident(ivar,2)
 np=-1
 if(vartyp eq 0) then np = npro
 if(vartyp eq 1) then np = 2
 if(vartyp eq 2) then np = 1
 if(vartyp eq 3) then np = 1
 if(vartyp eq 9) then np = 3
 if(vartyp eq 555) then np = 1


 ikeep = where(vconv ge a and vconv le b) 

 aa = fltarr(nconv,npro)
 aa(*)=0.0
 if(vartyp eq 0) then begin
  for igeom=0,ngeom-1 do begin
   for i=0,nconv-1 do begin
    aa(i,*)=aa(i,*)+kk(igeom*nconv+i,ioff:(ioff+npro-1))
   endfor
  endfor

  bb(*)=0.0
  ap = fltarr(npro)

  icut = where(lpref le lpress(0) or lpref ge lpress(npro-1),count)
  for i=0,nconv-1 do begin
   j = long((vconv(i)-vmin)/dv)
   ap(*)=aa(i,*)
   a1 = reverse(ap)
   b1 = interpol(a1,lpress,lpref)
   if(count gt 0) then b1(icut)=0.0
   bp = reverse(b1)
   bb(j,*)=bp(*)
  endfor

  print,'Enter title: '
  tname=''
  read,tname

  print,'Wavenumber(0) or wavelength(1)?'
  iwave=0
  read,iwave

  imagekkbar,bb,tname,iplotx,nplotx,iploty,nploty,$
   pmin,pmax,vmin,vmax,ips,ict,iwave

  iplotx = iplotx+1
  if(iplotx gt nplotx) then begin
   iploty = iploty+1
   iplotx = 1
  endif

 endif


 ioff = ioff+np

endfor

if(ips ne 0) then device,/close
 
end
