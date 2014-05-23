pro convertk,ntemp_in,nwave_in,temp_in,wave_in,k_in,$
             ntemp_out,nwave_out,temp_out,wave_out,k_out
; ***************************************************************
; Procedure to convert a k-table from one grid to another
;
; Pat Irwin	23/5/14
;
; ***************************************************************

k1=fltarr(ntemp_out,nwave_in)
tmp=fltarr(ntemp_in)
; Interpolate onto temperature array
for i=0,nwave_in-1 do begin
 tmp(*)=k_in(*,i)
 k1(*,i)=interpol(tmp,temp_in,temp_out)
endfor

; Check that for any temperatures outside the input ktable range, the k-values
; are set to the max/min temperature grid of the input array
ikeep = where(temp_out lt temp_in(0),nkeep)
if(nkeep gt 0) then for i=0,nkeep-1 do k1(ikeep(i),*)=k_in(0,*)
ikeep = where(temp_out gt temp_in(ntemp_in-1),nkeep)
if(nkeep gt 0) then for i=0,nkeep-1 do k1(ikeep(i),*)=k_in(ntemp_in-1,*)
 

;interpolate onto wavenumber array
k_out=fltarr(ntemp_out,nwave_out)
tmp=fltarr(nwave_in)
for i=0,ntemp_out-1 do begin
 tmp(*)=k1(i,*)
 k_out(i,*)=interpol(tmp,wave_in,wave_out)
endfor

; Check that for any wavenumbers outside the input ktable range, the k-values
; are set to zero.
ikeep = where(wave_out lt wave_in(0),nkeep)
if(nkeep gt 0) then for i=0,nkeep-1 do k_out(*,ikeep(i))=0.
ikeep = where(wave_out gt wave_in(nwave_in-1),nkeep)
if(nkeep gt 0) then for i=0,nkeep-1 do k_out(*,ikeep(i))=0.


return

end
