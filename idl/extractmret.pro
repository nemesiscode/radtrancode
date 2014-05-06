pro extractmret,filename,wkeep,spec,err,specf,irefl,pressA,heightA,idisoB,$
  xtemp,xcloud,xvmr,xradius,xtground
; ********************************************************************
; IDL program to read in and plot a Nemesis_3A .mre file
;
; Pat Irwin	12/2/04
; ********************************************************************

inpname = strcompress(filename + '.inp',/REMOVE_ALL)
ispace=0
openr,1,inpname
readf,1,ispace
close,1

retname = strcompress(filename + '.mre',/REMOVE_ALL)
apname = strcompress(filename + '.apr',/REMOVE_ALL)
prfname = strcompress(filename + '.prf',/REMOVE_ALL)
truname = strcompress(filename + '.rfa',/REMOVE_ALL)
refname = strcompress(filename + '.ref',/REMOVE_ALL)
aername = 'aerosol.prf'


; ******* Read in .prf file to get pressure grid *******
readprfhead4,prfname,npro,press,height,temp,molwt

; ******* Read in aerosol.prf file to get ncont *******
readaerprfhead,aername,ncont

; ******* Read in .apr file to get a priori measurement vector *******
readapriori,apname,npro,nvar,varident,varparam,nx,xa,erra


readprofilenew,truname,iformB,nplanetB,xlatB,nproB,ngasB,molwtB,idisoB,heightB,$
        pressB,tempB,vmrB
readprofilenew,prfname,iformA,nplanetA,xlatA,nproA,ngasA,molwtA,idisoA,heightA,$
        pressA,tempA,vmrA
readprofilenew,refname,iformC,nplanetC,xlatC,nproC,ngasC,molwtC,idisoC,heightC,$
        pressC,tempC,vmrC


xtemp=fltarr(5,npro)
xcloud=fltarr(4,3)
xradius=fltarr(4)
xtground=xradius
xradius(*)=-999.9
xtground(*)=-999.9
xvmr=fltarr(3,ngasB,npro)
xvmr(*)=-1.
for i=0,ngasB-1 do xvmr(0,i,*)=vmrB(*,i)
xtemp(4,*)=tempB(*)

openr,1,retname

 nspec = 1
 readf,1,nspec

 print,nspec
 iplot=1
 
 for ip = 0,iplot-1 do begin
   itmp = intarr(5)
   readf,1,itmp
   ngeom = itmp(1)
   nconv = itmp(2)
   nx = itmp(3)
   ny = itmp(4)
   latlon = fltarr(2)
   readf,1,latlon
   uname=''
   readf,1,uname
   head=''
   readf,1,head

   ydat = fltarr(7,ny)
   xdat = fltarr(4,nx)
   readf,1,ydat
   for i=1,2 do begin
    readf,1,head
   endfor
   istart=0
   for ivar=0,nvar-1 do begin
     print,ivar
     print,varident(ivar,*)
     for i=1,4 do begin
;     for i=1,44 do begin
      readf,1,head
     endfor
     itype = varident(ivar,2)
     print,'itype = ',itype
     case itype of
      0: np = npro
      1: np = 2
      2: np = 1
      3: np = 1
      4: np = 3
      6: np = 2
      8: np = 3
      9: np = 3
      10: np = 4
      11: np = 2
      12: np = 3
      13: np = 3
      14: np = 3
      15: np = 3
      555: np = 1
      888: np = varparam(ivar,0)
      889: np = 1
      999: np = 1
      777: np = 1
      666: np = 1
     endcase
     print,'np = ',np
     xdat1 = fltarr(6,np)
     readf,1,xdat1
     print,'istart = ',istart
     for j=0,np-1 do xdat(*,istart+j)=xdat1(2:5,j)

     istart=istart+np

   endfor

 endfor

close,1


if(irefl eq 1) then begin
  if(ispace eq 1) then begin
   ipfile = '/home/oxpln98/plan/irwin/radtrancode/trunk/raddata/sun_spec_echo.dat'
   print,'Sun file = ',ipfile
   openr,1,ipfile
   head=''
   for i=1,3 do readf,1,head
   solar = fltarr(2,600)
   readf,1,solar
   close,1
  endif
  AU=1.49598E13

  print,'Enter distance from Sun (AU) : '
  read,soldist

; Convert solar flux in units of uW cm-2 um-1
  area = 4*!pi*(soldist*AU)^2
  solar(1,*)=1e6*solar(1,*)/area

  print,'Enter zenith angle of Sun'
  read,zenkeep

endif

nconv1=ny
wkeep = fltarr(nconv1)
wkeep(*)=ydat(1,*)
spec = wkeep
err = wkeep
spec(*)=ydat(2,*)
err(*)=ydat(3,*)
specf=spec
specf(*)=ydat(5,*)

if(irefl eq 1) then begin
  solint = interpol(solar(1,*),solar(0,*),wkeep)
  sun = solint*cos(zenkeep*!pi/180.0)/!pi
  spec=spec/sun
  err=err/sun
  specf=specf/sun
endif


istart = 0
for ivar=0,nvar-1 do begin
  itype = varident(ivar,2)
  print,'ivar, itype',ivar,itype
  print,'varident = ',varident(ivar,*)

  case itype of
      0: np = npro
      1: np = 2
      2: np = 1
      3: np = 1
      4: np = 3
      6: np = 2
      8: np = 3
      9: np = 3
      12: np = 3
      13: np = 3
      14: np = 3
      15: np = 3
      555: np = 1
      888: np = varparam(ivar,0)
      889: np = 1
      999: np = 1
      777: np = 1
      666: np = 1
  endcase
  jgas=-1

  if(itype eq 0) then begin
     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))
 
     if(varident(ivar,0) eq 0) then begin
      xtemp(0,*)=xa(*)
      xtemp(1,*)=erra(*)
      xtemp(2,*)=xn(*)
      xtemp(3,*)=errn(*)
     endif

  endif else begin
   if(itype eq 9) then begin

     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))
     xcloud(0,*)=xa(*)
     xcloud(1,*)=erra(*)
     xcloud(2,*)=xn(*)
     xcloud(3,*)=errn(*)

   endif else begin

     if(itype eq 3) then begin
      jgas=-1
      for igas=0,ngasB-1 do begin
       if(varident(ivar,0) eq idisoB(0,igas) and varident(ivar,1) eq idisoB(1,igas))then begin
        jgas=igas
       endif
      endfor

      if(jgas ge 0) then begin
       xx1=xdat(2,istart)
       xx2=xdat(3,istart)
       xvmr(1,jgas,*)=xx1*vmrC(*,jgas)
       xvmr(2,jgas,*)=xx2*vmrC(*,jgas)
      endif

     endif
   endelse

  endelse

  if(itype eq 555) then begin
   xradius(*)=xdat(*,istart)
  endif
  if(itype eq 999) then begin
   xtground(*)=xdat(*,istart)
  endif

  istart = istart+np

endfor

return

end
