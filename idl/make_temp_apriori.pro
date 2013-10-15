;#######################################
;# idl routine to generate a continuous temperature profile
;#
;# Original by: bejamin mort 	2/02/04
;#	
;# Updated by:  Pat Irwin	28/7/04
;# Revised:	Pat Irwin	30/7/04
;#######################################

;read in pressures from .ref file
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Print,'Enter name of reference .ref file : '
filename=''
read,filename
;filename = 'cirsret.ref'
openr,1,filename

comment =''
head1: readf,1,comment
       print,comment
       if(strmid(comment,0,1) eq '#') then goto,head1

iform=1
reads,comment,iform

print,comment,'iform = ',iform

tmp=fltarr(5)
if(iform eq 1) then tmp = fltarr(4)
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

print,'Enter range of pressures to plot over : '
pmax = 1.
pmin = 1.
read,pmax,pmin

!p.multi=0
plot_io,temp,press,yrange=[pmax,pmin]


print,'Enter temperature error at top of atmosphere'
err1 = 1.0
read,err1
print,'Enter pressure level of unit optical depth'
pod1 = 1.0
read,pod1

print,'Enter CLEN'
read,CLEN


print,'Enter any temperature offset : '
toff=1.
read,toff

;calculate profile
;~~~~~~~~~~~~~~~~~
profile = fltarr (3,np)
profile(0,*) = press(*)
profile(1,*) = temp(*)+toff
odepth = press/pod1
profile(2,*)=err1*exp(-odepth)

plot_io,profile(1,*),press,yrange=[pmax,pmin]
oplot,profile(1,*)+profile(2,*),press,linestyle=1
oplot,profile(1,*)-profile(2,*),press,linestyle=1


;write temp profile
;~~~~~~~~~~~~~~~~~

print,'enter name of output file (.dat)'
read,filename
;print,'--------------------------------'
openw,2,filename
printf,2,np,CLEN
printf,2,profile
close,2


end


