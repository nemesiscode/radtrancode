;#######################################
;# idl routine to generate a continuous parah2 profile
;#
;# Original by: bejamin mort 	2/02/04
;#	
;# Updated by:  Pat Irwin	28/7/04
;# Revised:	Pat Irwin	30/7/04
;#######################################

;read in pressures from .ref file
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Print,'Enter name of reference parah2 .ref file : '
filename=''
read,filename
;filename = 'cirsret.ref'
openr,1,filename
 head=''
 readf,1,head
 np = 1
 readf,1,np
 parah2 = fltarr(2,np)
 readf,1,parah2
close,1


print,'Enter range of pressures to plot over : '
pmax = 1.
pmin = 1.
read,pmax,pmin

!p.multi=0
plot_io,parah2(1,*),press,yrange=[pmax,pmin]

print,'Enter paraH2 error'
ferr = 1.0
read,ferr

print,'Enter CLEN'
read,CLEN


;calculate profile
;~~~~~~~~~~~~~~~~~
profile = fltarr (3,np)
profile(0,*) = press(*)
profile(1,*) = parah2(1,*)
profile(2,*) = ferr

oplot,profile(1,*)+profile(2,*),press,linestyle=1
oplot,profile(1,*)-profile(2,*),press,linestyle=1

ikeep = where(profile(1,*) lt 1e-30,count)
if(count gt 0) then profile(1,ikeep)=1e-30
ikeep = where(profile(2,*) lt 1e-36,count)
if(count gt 0) then profile(2,ikeep)=1e-36

;write parah2 profile
;~~~~~~~~~~~~~~~~~

print,'enter name of output file (.dat)'
read,filename
;print,'--------------------------------'
openw,2,filename
printf,2,np,CLEN
printf,2,profile
close,2


end


