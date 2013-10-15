;#######################################
;# idl routine to generate a continuous and NH3 or PH3 profile
;#
;# Original by: bejamin mort 	2/02/04
;#	
;# Updated by:  Pat Irwin	28/7/04
;# Revised:	Pat Irwin	30/7/04
;#######################################

;--- NH3 ----
pknee = 0.7
deep_vmr = 2.19e-4
deep_vmr_err = 1e-4
fsh= .15
fsh_err = .1
CLEN = 1.5

;--- PH3 ----
;pknee = 1
;deep_vmr = 6e-7
;deep_vmr_err = 3e-7
;fsh= .3
;fsh_err = .2
;CLEN = 1.5

;read in pressures from .ref file
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Print,'Enter name of reference .ref file : '
filename=''
read,filename
;filename = 'cirsret.ref'
openr,1,filename

comment =''
head1: readf,1,comment
       if(strmid(comment,0,1) eq '#') then goto,head1

tmp = fltarr(5)
readf,1,tmp
np = long(tmp(2))
ngas = long(tmp(3))
print,tmp
idiso = lonarr(2,ngas)
readf,1,idiso
print,idiso

height = fltarr(np)
press = height
temp = height
vmr = fltarr(ngas,np)
remgas = ngas
ioff = 0
ncond = 3
loop1: ncol = min([6,ncond + remgas])
       data = fltarr(ncol,np)
       readf,1,comment
       readf,1,data
 
       if(ncond eq 3) then begin
        height(*)=data(0,*)
        press(*)=data(1,*)
        temp(*)=data(2,*)
        for igas = 1,ncol-3 do vmr(ioff+igas-1,*)=data(3+igas-1,*)
       endif else for igas = 1,ncol do vmr(ioff+igas-1,*)=data(igas-1,*)
       ioff = ioff + ncol-ncond
       remgas = remgas - (ncol-ncond)
       print,ioff,remgas

       ncond = 0
       if(remgas gt 0) then goto,loop1
close,1

print,'Which gas profile do you want to select : '
for i=0,ngas-1 do print,i,idiso(*,i)
read,igas

print,'Enter range of pressures to plot over : '
pmax = 1.
pmin = 1.
read,pmax,pmin

!p.multi=0
plot_oo,vmr(igas,*),press,yrange=[pmax,pmin]


print,'Enter pknee'
read,pknee
print,'Enter deep_vmr,deep_vmr_err'
read,deep_vmr,deep_vmr_err
print,'Enter fsh,fsh_err'
read,fsh,fsh_err
print,'Enter CLEN'
read,CLEN


print,'----------------------------'
print, 'pknee -- deep vmr -- fsh -- correlation length'
print,pknee,deep_vmr,fsh,CLEN
print, '----------------------------'



;calculate profile
;~~~~~~~~~~~~~~~~~
profile = fltarr (3,np)
profile(0,*) = press(*)
for i=0,np-1 do begin
 if (profile(0,i) ge pknee) then begin
  profile(1,i) = deep_vmr
  profile(2,i) = deep_vmr_err
 endif

 if (profile(0,i) lt pknee) then begin
  p_po=profile(0,i)/pknee
  profile(1,i) = deep_vmr * (p_po)^((1-fsh)/fsh)
  profile(2,i) = p_po^((1-fsh)/fsh)*deep_vmr_err+((-p_po^(1/fsh)*alog(p_po))/(fsh^2*p_po))*deep_vmr*fsh_err
 endif

endfor 


xvmr = press
xvmr_error = press
xvmr(*)=profile(1,*)
xvmr_error(*) = profile(2,*)

ikeep = where(press ge pmin and press le pmax)

vmr1 = exp(alog(xvmr)-xvmr_error/xvmr)
vmr2 = exp(alog(xvmr)+xvmr_error/xvmr)

a = min(vmr1(ikeep))
b = max(vmr2(ikeep))

oplot,xvmr,press

oplot,vmr1,press,linestyle=1
oplot,vmr2,press,linestyle=1

ikeep = where(profile(1,*) lt 1e-36,count)
if(count gt 0) then profile(1,ikeep)=1e-36
if(count gt 0) then profile(2,ikeep)=1e-30
;ikeep = where(profile(2,*) lt 1e-30,count)
;if(count gt 0) then profile(2,ikeep)=1e-30

;write gas profile
;~~~~~~~~~~~~~~~~~

print,'enter name of output file (.dat)'
read,filename
;print,'--------------------------------'
openw,2,filename
printf,2,np,CLEN
printf,2,profile
close,2


end


