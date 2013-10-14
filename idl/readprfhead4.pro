pro readprfhead4,prfname,npro,press,height,temp,molwt
; ******************************************************************
; IDL routine to read in a .prf file to get the pressure grid
; Input variables
;    prfname	character	filename
;
; Output variables
;    npro	integer		Number of levels in .prf file
;    press(npro) real		Pressure levels
;    height(npro) real		Height levels
;    temp(npro)  real		Temperatures
;    molwt	real		Molecular weight of atmosphere (kg/mol)
;
; Pat Irwin	28/10/03
; ******************************************************************

openr,4,prfname
iform=1
head=''
readf,4,iform
tmp=fltarr(4)
if(iform eq 0) then tmp = fltarr(5)
readf,4,tmp
npro = long(tmp(2))
ngas = long(tmp(3))
molwt=-1.
if(iform eq 0) then molwt=tmp(4)*1e-3
for i=1,ngas+1 do readf,4,head
ncol = 3+ngas
if(ncol gt 6) then ncol=6
tdat = fltarr(ncol,npro)
readf,4,tdat
close,4

press = fltarr(npro)
height=press
; convert atm to bar for pressure
press(*)=tdat(1,*)*1.013
temp=fltarr(npro)
temp(*)=tdat(2,*)
height(*)=tdat(0,*)

return

end

