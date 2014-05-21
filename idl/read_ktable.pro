; *********************************************************************
; Procedure to read in a Nemesis k-table into IDL and examine the contents
;
; Pat Irwin  21/5/14
;
; *********************************************************************

filename=''
print,'Enter filename : '
read,filename
openr,1,filename,/rawio
irec0=long(1)
npoint=irec0
vmin=1.0
delv=1.0
fwhm=1.0
np=irec0
nt=irec0
ng=irec0
idgas=irec0
isogas=irec0

readu,1,irec0
readu,1,npoint
readu,1,vmin
readu,1,delv
readu,1,fwhm
readu,1,np
readu,1,nt
readu,1,ng
readu,1,idgas
readu,1,isogas

print,'irec0,npoint,vmin,delv,fwhm = ',irec0,npoint,vmin,delv,fwhm
print,'idgas,isogas = ', idgas,isogas

x=1.
gord=fltarr(ng)
delg=gord
print,' '
print,'ng = ',ng
for i=0,ng-1 do begin
 readu,1,x
 gord(i)=x
endfor 
for i=0,ng-1 do begin
 readu,1,x
 delg(i)=x 
endfor

for i=0,ng-1 do begin
 print,i+1,gord(i),delg(i)
endfor

; skip over two records
for i=0,1 do readu,1,x

press=fltarr(np)
temp=fltarr(nt)
print,' '
print,'np = ',np
for i=0,np-1 do begin
 readu,1,x
 press(i)=x
 print,i+1,press(i)
endfor 
print,' '
print,'nt = ',nt
for i=0,nt-1 do begin
 readu,1,x
 temp(i)=x 
 print,i+1,temp(i)
endfor

; Number of records read so far
nrec = 10+2*ng+2+np+nt


; Read in wavelengths if non-uniform grid
vwave=fltarr(npoint)
if(delv lt 0.0) then begin
 for i=0,npoint-1 do begin
  readu,1,x
  vwave(i)=x
 endfor
 nrec=nrec+npoint
endif else vwave=vmin+delv*findgen(npoint)

;print,' '
;print,'vwave : '
;for i=0,npoint-1 do print,i+1,vwave(i)

; Number of records to skip to get to first record of ktable
jrec = irec0-nrec-1
if(jrec gt 0) then for i=1,jrec do readu,1,x

ktable=fltarr(npoint,np,nt,ng)

for ipoint=0,npoint-1 do begin
 for j=0,np-1 do for k=0,nt-1 do for i=0,ng-1 do begin
  readu,1,x
  ktable(ipoint,j,k,i)=x
 endfor
endfor
close,1

loop:
print,'temperature(1-nt) and g-ordinate(1-ng) to examine : '
read,it,ig

if(it lt 1 or it gt nt) then goto,loop 
if(ig lt 1 or ig gt ng) then goto,loop 

tmp=fltarr(npoint,np)
for i=0,np-1 do tmp(*,i)=ktable(*,i,it-1,ig-1)

surface,tmp,vwave,alog10(press)

print,'Plot another combination (Y/N):'
ans=''
read,ans
if(strupcase(ans) eq 'Y') then goto,loop

end
