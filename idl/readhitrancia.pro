pro readhitrancia,filename,ntemp,nwave,temp,wave,data

openr,1,filename
head=''
nwave=1
v1=1.
v2=1.
x1=1.
for itemp=0,ntemp-1 do begin
 readf,1,head
; print,head
 x2 = strmid(head,20,99-20)
 reads,x2,v1,v2,nwave,x1
 if(itemp eq 0)then begin
  tmp=dblarr(2,nwave)
  wave=fltarr(nwave)
  data=fltarr(ntemp,nwave)
  temp=fltarr(ntemp)
 endif
 temp(itemp)=x1
 readf,1,tmp
; Hitran format CIA is in units of cm5 mol-2. Convert to cm-1 amagat-2 using
; conversion factor defined by Richard et al. 2012.
 data(itemp,*)=float(tmp(1,*)/1.385e-39)
 if(itemp eq 0) then wave(*)=tmp(0,*)
endfor
close,1

return

end

