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
; print,ivar+1,itmp
 varident(ivar,*)=itmp(*)
 itype = itmp(2)

; print,'ITYPE = ',ITYPE
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


 if(itype eq 29) then begin
  apfile=''
  readf,5,apfile
  apfile = strcompress(apfile,/REMOVE_ALL)
  openr,6,apfile
   npro1=1
   nloc=1
   readf,6,nloc,npro1,clen
   if(npro1 ne npro)then begin
    print,'Error in readapriori.pro'
    stop
   endif
   np = nloc*npro
   data = fltarr(3,np)
   readf,6,data
  close,6
  xa(istart:istart+np-1)=data(1,*)
  erra(istart:istart+np-1)=data(2,*)
 endif 

 if(itype eq 1 or itype eq 6 or itype eq 7) then begin
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

 if(itype eq 30) then begin
  apfile=''
  readf,5,apfile
  apfile = strcompress(apfile,/REMOVE_ALL)
  openr,6,apfile
   nlong=1
   nlevel=1
   clen1=1.
   clen2=2.
   head=''
   readf,6,head
   if(strmid(head,0,2) eq 'Ex') then begin
    readf,6,nlong,nlevel,clen1,clen2
   endif else reads,head,nlong,nlevel,clen1,clen2

   varparam(ivar,0)=nlong*nlevel
   varparam(ivar,1)=nlevel
   data = fltarr(nlong+2,nlevel)
   readf,6,data
   xpc =1.0
   readf,6,xpc
  close,6
  varparam(ivar,2:2+nlevel-1)=data(0,*)
  varparam(ivar,2+nlevel)=xpc

  np = nlong*nlevel

  for ilevel=0,nlevel-1 do begin
   for ilong=0,nlong-1 do begin
    k = ilong*nlevel+ilevel
    xa(istart+k)=data(ilong+1,ilevel)
    erra(istart+k)=data(nlong+1,ilevel)
   endfor
  endfor
 endif


 if(itype eq 31) then begin
  apfile=''
  readf,5,apfile
  apfile = strcompress(apfile,/REMOVE_ALL)
  openr,6,apfile
   nlong=1
   clen2=2.
   readf,6,nlong,clen2
   varparam(ivar,0)=nlong
   data = fltarr(nlong+1)
   readf,6,data
   xpc=1.0
   readf,6,xpc
   varparam(ivar,1)=xpc
  close,6
  
  np = nlong

  for k=0,nlong-1 do begin
    xa(istart+k)=data(k)
    erra(istart+k)=data(nlong)
  endfor

 endif

 if(itype eq 2 or itype eq 3) then begin
  np = 1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
 endif 

 if(itype eq 8 or itype eq 9 or itype eq 32 or itype eq 45) then begin
  np = 3
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+2)=tmp(0,*)
  erra(istart:istart+2)=tmp(1,*)
 endif 
 
 if(itype eq 37) then begin
  np = 1
  head=''
  readf,5,head
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
 endif

 if(itype eq 42) then begin
  np = 3
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+2)=tmp(0,*)
  erra(istart:istart+2)=tmp(1,*)
  head=''
  readf,5,head
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

 if(itype eq 21) then begin
  np = 2
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+1)=tmp(0,*)
  erra(istart:istart+1)=tmp(1,*)
  head=''
  readf,5,head  
 endif 

 if(itype eq 40) then begin
  np = 2
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+1)=tmp(0,*)
  erra(istart:istart+1)=tmp(1,*)
  tmp=fltarr(3)
  readf,5,tmp
  varparam(ivar,0:2) = tmp(*)
 endif 

 if(itype ge 12 and itype le 15) then begin
  np = 3
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+2)=tmp(0,*)
  erra(istart:istart+2)=tmp(1,*)
 endif 

 if(itype eq 19) then begin
  np = 4
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+3)=tmp(0,*)
  erra(istart:istart+3)=tmp(1,*)
 endif 

 if(itype eq 39) then begin
  np = 1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
  readf,5,tmp
  varparam(ivar,0:1) = tmp(*)
 endif 

 if(itype eq 22) then begin
  np = 5
  tmp = fltarr(2,np)
  tmp1=fltarr(2)
  for i=0,4 do begin
   readf,5,tmp1
   tmp(*,i)=tmp1(*)
  endfor
  xa(istart:istart+4)=tmp(0,*)
  erra(istart:istart+4)=tmp(1,*)
 endif 

 if(itype eq 555 or itype eq 333) then begin
  np = 1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart)=tmp(0)
  erra(istart)=tmp(1)
 endif 

 
 if(itype eq 25)then begin
  x=1.
  np=1
  readf,5,np,x
  varparam(ivar,0) = np
  tmp = fltarr(2,np)
  readf,5,tmp
  xa(istart:istart+np-1)=tmp(0,*)
  erra(istart:istart+np-1)=tmp(1,*)
 endif

 if(itype eq 889)then begin
  np = 1
  tmp = fltarr(2)
  readf,5,tmp
  xa(istart:istart+np-1)=tmp(0)
  erra(istart:istart+np-1)=tmp(1)
 endif

 if(itype eq 888)then begin
  readf,5,np
  varparam(ivar,0) = np
  tmp = fltarr(3,np)
  readf,5,tmp
  xa(istart:istart+np-1)=tmp(1,*)
  erra(istart:istart+np-1)=tmp(2,*)
 endif

 if(itype eq 887)then begin
  readf,5,np
  varparam(ivar,0) = np
  tmp = fltarr(3,np)
  readf,5,tmp
  xa(istart:istart+np-1)=tmp(1,*)
  erra(istart:istart+np-1)=tmp(2,*)
 endif

 if(itype eq 222 or itype eq 226)then begin
  np=8
  tmp = fltarr(2,np)
  xx=fltarr(2)
  for k=0,7 do begin
    readf,5,xx
    tmp(*,k)=xx(*)
  endfor
  if(itype eq 222) then begin
   xtmp=fltarr(5)
   readf,5,xtmp
   varparam(ivar,0:4) = xtmp(*)
  endif

  xa(istart:istart+np-1)=tmp(0,*)
  erra(istart:istart+np-1)=tmp(1,*)

 endif

 if(itype eq 225)then begin
  np=11
  tmp = fltarr(2,np)
  xx=fltarr(2)
  for k=0,10 do begin
    readf,5,xx
    tmp(*,k)=xx(*)
  endfor
  xtmp=fltarr(5)
  readf,5,xtmp
  varparam(ivar,0:4) = xtmp(*)

  xa(istart:istart+np-1)=tmp(0,*)
  erra(istart:istart+np-1)=tmp(1,*)

 endif


 if(itype eq 223 or itype eq 224)then begin
  np=9
  tmp = fltarr(2,np)
  xx=fltarr(2)
  for k=0,8 do begin
    readf,5,xx
    tmp(*,k)=xx(*)
  endfor
  xtmp=fltarr(5)
  readf,5,xtmp
  varparam(ivar,0:4) = xtmp(*)

  xa(istart:istart+np-1)=tmp(0,*)
  erra(istart:istart+np-1)=tmp(1,*)

 endif

 if(itype eq 444)then begin
  cloudfil=''
  readf,5,cloudfil
  openr,10,cloudfil
   np=2
   xx=fltarr(2)
   tmp=fltarr(2,np)
   for k=0,1 do begin
    readf,10,xx
    tmp(*,k)=xx(*)
   endfor
;   print,tmp
   xa(istart:istart+1)=tmp(0,*)
   erra(istart:istart+1)=tmp(1,*)
   tmp=fltarr(2)
   readf,10,tmp
   nlam=tmp(0)
   clen=tmp(1)
;   print,'nlam = ',nlam
   varparam(ivar,0:1) = tmp(*)
   readf,10,tmp
   varparam(ivar,2:3) = tmp(*)
   xlam=1.
   readf,10,xlam
   varparam(ivar,4) = xlam
   tmp=fltarr(3,nlam)
   readf,10,tmp
   if(clen gt 0) then begin
    xa(istart+2:istart+2+nlam-1)=tmp(1,*)
    erra(istart+2:istart+2+nlam-1)=tmp(2,*)
    np=np+nlam
   endif else begin
    xa(istart+2)=tmp(1,0)
    erra(istart+2)=tmp(2,0)
    np=np+1   
   endelse
   walb=fltarr(nlam)
   walb(*)=tmp(0,*)
  close,10
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
