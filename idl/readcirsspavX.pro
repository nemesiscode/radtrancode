pro readcirsspavX,filename,fwhm,xlat,xlon,ngeom,nav,nconv,wave,$
       angles,spec,error
; **********************************************************************
; IDL procedure to read in a Nemesis .spc file
; Input variables
;    filename	character	filename
;
; Output variables
;    fwhm	real		FWHM of spectrum
;    xlat	real		latitude
;    xlon	real		longitude
;    ngeom	integer		Number of observations
;    nav(ngeom)	integer		FOV averaging 
;    nconv	integer		Length of each sub-spectrum
;    wave(nconv) real		Measured wavenumbers
;    angles(6,ngeom,20) real	Observation angles
;    spec(ngeom,nconv) real	Measured spectra
;    error(ngeom,nconv) real	Measured errors
;
; Pat Irwin	29/10/03
; **********************************************************************

openr,1,filename
tmp = fltarr(4)
readf,1,tmp
fwhm=tmp(0)
xlat=tmp(1)
xlon=tmp(2)
ngeom = long(tmp(3))
nav = lonarr(ngeom)
angles = fltarr(6,ngeom,20)

for igeom=0,ngeom-1 do begin
 readf,1,nconv
 if(igeom eq 0) then begin
  spe = fltarr(3,nconv)
  wave = fltarr(nconv)
  spec = fltarr(ngeom,nconv)
  error = fltarr(ngeom,nconv)
  ang1 = fltarr(6)
 endif
 i1=1
 readf,1,i1
 nav(igeom)=i1
 for iav=0,i1-1 do begin
  readf,1,ang1
  angles(*,igeom,iav)=ang1(*)
 endfor
 readf,1,spe
 if(igeom eq 0) then wave = spe(0,*)
 spec(igeom,*)=spe(1,*)
 error(igeom,*)=spe(2,*)
endfor

close,1

return

end
