pro readapriori,apname,npro,nvar,varident,varparam,nx,xa,erra,walb
; ********************************************************************
; Procedure to read in a Nemesis .apr file
; Input variables
;   apname	character	filename
;   npro	integer		Number of vertical levels
;
; Output variables
;   nvar	integer		Number of variable species
;   varident(mvar,3) integer    Identity of each variable species
;   varparam(mvar,nparam) real  Extra variables as needed
;   nx		integer		Length of measurement vector
;   xa(mx)	real		Measurement a priori vector
;   err(mx)	real		A priori errors
;
; Pat Irwin	29/10/03
; ********************************************************************

mx = 400
mvar = 12
np = 1
mparam = 405
xa = fltarr(mx)
erra = fltarr(mx)
varident = intarr(mvar,3)
varparam = fltarr(mvar,mparam)
head=''

;print,mvar,size(varident)

openr,5,apname
head=' '
readf,5,head
nvar=1
readf,5,nvar
istart=0

for ivar=0,nvar-1 do begin
 itmp = intarr(3)
 readf,5,itmp
 varident(ivar,*)=itmp(*)
 itype = itmp(2)

 print,'ITYPE = ',ITYPE
 if(itype eq 0) then begin
  np = npro
  apfile=''
  readf,5,apfile
  apfile = strcompress(apfile,/REMOVE_ALL)
  openr,6,apfile
   npro1=1
   readf,6,npro1,clen
   if(npro1 ne npro)then begin
    print,'Error in readapriori.pro'
    stop
   endif
   data = fltarr(3,npro)
   readf,6,data
  close,6
  xa(istart:istart+np-1)=data(1,*)
  erra(istart:istart+np-1)=data(2,*)
 endif 

 if(itype eq 1 or itype eq 6) then begin
  np = 2
  readf,5,xp
  varparam(ivar,0)=xp
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
  readf,5,tmp
  xa(istart+1)=tmp(0)
  erra(istart+1)=tmp(1)
 endif 

 if(itype eq 2 or itype eq 3) then begin
  np = 1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
 endif 

 if(itype eq 8 or itype eq 9) then begin
  np = 3
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+2)=tmp(0,*)
  erra(istart:istart+2)=tmp(1,*)
 endif 

 if(itype eq 10) then begin
  np = 4
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+3)=tmp(0,*)
  erra(istart:istart+3)=tmp(1,*)
  head=''
  readf,5,head
 endif 

 if(itype eq 11) then begin
  np = 2
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+1)=tmp(0,*)
  erra(istart:istart+1)=tmp(1,*)
  head=''
  readf,5,head  
 endif 

 if(itype ge 12 and itype le 15) then begin
  np = 3
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+2)=tmp(0,*)
  erra(istart:istart+2)=tmp(1,*)
 endif 

 if(itype eq 555) then begin
  np = 1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
 endif 

 if(itype eq 888)then begin
  readf,5,np
  varparam(ivar,0) = np
  tmp = fltarr(3,np)
  readf,5,tmp
  xa(istart:istart+np-1)=tmp(1,*)
  erra(istart:istart+np-1)=tmp(2,*)
 endif

 if(itype eq 444)then begin
  cloudfil=''
  readf,5,cloudfil
  openr,10,cloudfil
   np=2
   tmp=fltarr(2,np)
   readf,10,tmp
   print,tmp
   xa(istart:istart+1)=tmp(0,*)
   erra(istart:istart+1)=tmp(1,*)
   tmp=fltarr(2)
   readf,10,tmp
   nlam=tmp(0)
   print,'nlam = ',nlam
   varparam(ivar,0:1) = tmp(*)
   readf,10,tmp
   varparam(ivar,2:3) = tmp(*)
   xlam=1.
   readf,10,xlam
   varparam(ivar,4) = xlam
   tmp=fltarr(3,nlam)
   readf,10,tmp
   xa(istart+2:istart+2+nlam-1)=tmp(1,*)
   erra(istart+2:istart+2+nlam-1)=tmp(2,*)
   walb=fltarr(nlam)
   walb(*)=tmp(0,*)
  close,10
  np=np+nlam
 endif

 if(itype eq 666)then begin
  np=1
  x1=1.0
  readf,5,x1
  varparam(ivar,0) = x1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
 endif

 
 if(itype eq 999 or itype eq 777) then begin
  np = 1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
 endif 

 istart = istart+np

endfor

close,5

nx = istart

return 

end
