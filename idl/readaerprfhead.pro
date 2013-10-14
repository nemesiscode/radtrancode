pro readaerprfhead,prfname,ncont
; ******************************************************************
; IDL routine to read in a aerosol.prf file to get the number of 
; cloud types
; Input variables
;    prfname	character	filename
;
; Output variables
;    ncont	integer		Number of cloud types
;
; Pat Irwin	12/2/04
; ******************************************************************

openr,1,prfname
 head=''
 readf,1,head
 tmp = intarr(2)
 readf,1,tmp
 ncont = tmp(1)
close,1

return

end

