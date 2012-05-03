; ********************************************************************
; IDL program to read in and plot a Nemesis_3A .mre file
;
; Pat Irwin	12/2/04
; ********************************************************************

!p.position=0

print,'Enter Nemesis run root name'
filename=' '
read,filename

print,'Plot to window(0), ps(1) or eps(2)?'
ips=0
read,ips

retname = strcompress(filename + '.mre',/REMOVE_ALL)
apname = strcompress(filename + '.apr',/REMOVE_ALL)
prfname = strcompress(filename + '.prf',/REMOVE_ALL)
aername = 'aerosol.prf'

; ******* Read in .prf file to get pressure grid *******
readprfhead4,prfname,npro,press,height,temp,molwt

; ******* Read in aerosol.prf file to get ncont *******
readaerprfhead,aername,ncont

; ******* Read in .apr file to get a priori measurement vector *******
readapriori,apname,npro,nvar,varident,varparam,nx,xa,erra

openr,1,retname

 nspec = 1
 readf,1,nspec

 print,nspec

 print,'There are ',nspec,' retrievals in total'
 print,'Enter the number you want to inspect (1-nspec)'
 read,iplot
 
 for ip = 0,iplot-1 do begin
   itmp = intarr(5)
   readf,1,itmp
   ngeom = itmp(1)
   nconv = itmp(2)
;   print,'nconv = ',nconv
;   print,'OK? (Y/N)'
;   ans=''
;   read,ans
;   if(ans ne 'y' and ans ne 'Y')then begin
;    print,'Enter new value : '
;    read,nconv
;   endif
   nx = itmp(3)
   ny = itmp(4)
   latlon = fltarr(2)
   readf,1,latlon
   head=''
   for i=1,2 do begin
    readf,1,head
   endfor

   ydat = fltarr(7,ny)
   xdat = fltarr(4,nx)
   readf,1,ydat
   for i=1,2 do begin
    readf,1,head
   endfor
   istart=0
   for ivar=0,nvar-1 do begin
     for i=1,4 do begin
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
      888: np = varparam(ivar,0)
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

print,'Wavenumber range = ',min(ydat(1,*)),max(ydat(1,*))
print,'wavenumber (0) or wavelength (1)?'
iwav=0
read,iwav

print,'Enter required plot range : '
xr = fltarr(2)
read,xr

print,'Plot radiance (0) or reflectivity (1) : '
irefl=0
read,irefl

 if(irefl eq 1) then begin
  if(iwav eq 1) then begin 
   openr,3,'/home/oxpln98/plan/irwin/baldrick/raddata/solar.dat'
   for i=1,3 do readf,3,head
  endif else begin
   openr,3,'/home/oxpln98/plan/irwin/baldrick/raddata/solWNum.dat'
   for i=1,5 do readf,3,head
  endelse
  head=' '
  solar=fltarr(2,78)
  readf,3,solar
  close,3

  if(iwav eq 1) then begin
;  input units are W m-2 nm-1 at 1AU. Convert to uW cm-2 um-1'
   solar(1,*)=solar(1,*)*1e-4*1e3*1e6
  endif else begin
;  input units are W m-2 (cm-1)-1 at 1AU. Convert to nW cm-2 (cm-1)-1'
   solar(1,*)=solar(1,*)*1e-4*1e9   

   tmp=fltarr(78)
   tmp(*)=solar(0,*)
   solar(0,*)=reverse(tmp)
   tmp(*)=solar(1,*)
   solar(1,*)=reverse(tmp)
  endelse
 
  print,'Enter distance from Sun (AU) : '
  read,soldist
;  soldist=19.2

;  print,'Enter solar zenith angle : '
;  read,zenkeep
   zenkeep=0.
 endif

if(ips eq 0) then begin
  set_plot,'x'
  device,retain=2
  device,decomposed=0
  window,0,title='Spectrum',xsize=600,ysize=600
endif else begin
 set_plot,'ps'
 if(ips eq 1) then begin
   device,filename='spectrum.ps',encapsulated=0 
 endif else device,filename='spectrum.eps',/encapsulated 
endelse

nconv1=ny/ngeom
!p.multi=[0,1,2]
i1=0
i2=nconv1-1
wkeep = fltarr(nconv1)
wkeep(*)=ydat(1,i1:i2)
spec = wkeep
err = wkeep
spec(*)=ydat(2,i1:i2)
err(*)=ydat(3,i1:i2)
specf=spec
specf(*)=ydat(5,i1:i2)

if(irefl eq 1) then begin
 solint = interpol(solar(1,*),solar(0,*),wkeep)
 sun = solint/(!pi*soldist^2)
 spec=spec/sun
 err=err/sun
 specf=specf/sun
 yname = 'Reflectivity'
endif else begin
 if(iwav eq 0) then begin
  yname = 'Radiance (nW cm!E-2!N sr!E-1!N (cm!E-1!N)!E-1!N)'
 endif else yname = 'Radiance (!4l!xW cm!E-2!N sr!E-1!N !4l!xm!E-1!N)'
endelse

y2=max([spec+err,specf])
y1=0.

print,'y1,y2 = ',y1,y2

xname = 'Wavenumbers (cm!E-1!N)'
if(iwav eq 1) then xname = 'Wavelength (!4l!xm)'

plot,wkeep,spec,xtitle=xname,$
ytitle=yname,xrange=xr,yrange=[y1,y2],xstyle=1
oplot,wkeep,spec+err,linestyle=1
oplot,wkeep,spec-err,linestyle=1
oplot,wkeep,specf,linestyle=2

for igeom=1,ngeom-1 do begin
 i1 = igeom*nconv1
 i2 = i1+nconv1-1

 wkeep(*)=ydat(1,i1:i2)
 spec(*)=ydat(2,i1:i2)
 err(*)=ydat(3,i1:i2)
 specf(*)=ydat(5,i1:i2)

 if(irefl eq 1) then begin
  solint = interpol(solar(1,*),solar(0,*),wkeep)
  sun = solint*cos(zenkeep*!pi/180.0)/(!pi*soldist^2)
  spec=spec/sun
  err=err/sun
  specf=specf/sun
 endif


 oplot,wkeep,spec
 oplot,wkeep,spec+err,linestyle=1
 oplot,wkeep,spec-err,linestyle=1
 oplot,wkeep,specf,linestyle=2

endfor

i1=0
i2=nconv1-1

wkeep(*)=ydat(1,i1:i2)
spec(*)=ydat(2,i1:i2)
err(*)=ydat(3,i1:i2)
specf(*)=ydat(5,i1:i2)

if(irefl eq 1) then begin
 solint = interpol(solar(1,*),solar(0,*),wkeep)
 sun = solint*cos(zenkeep*!pi/180.0)/(!pi*soldist^2)
 spec=spec/sun
 err=err/sun
 specf=specf/sun
 yname = 'Reflectivity'
endif else yname = 'Radiance (nW cm!E-2!N sr!E-1!N (cm!E-1!N)!E-1!N)'

y2=max([specf-spec,err])
y1=min([specf-spec,-err])


plot,wkeep,specf-spec,xtitle='Wavenumbers (cm!E-1!N)',$
ytitle=yname,xrange=xr,yrange=[y1,y2],xstyle=1
oplot,wkeep,err,linestyle=1
oplot,wkeep,-err,linestyle=1

for igeom=0,ngeom-1 do begin
 i1 = igeom*nconv1
 i2 = i1+nconv1-1

 wkeep(*)=ydat(1,i1:i2)
 spec(*)=ydat(2,i1:i2)
 err(*)=ydat(3,i1:i2)
 specf(*)=ydat(5,i1:i2)

 if(irefl eq 1) then begin
  solint = interpol(solar(1,*),solar(0,*),wkeep)
  sun = solint*cos(zenkeep*!pi/180.0)/(!pi*soldist^2)
  spec=spec/sun
  err=err/sun
  specf=specf/sun
 endif

 oplot,wkeep,specf-spec
 oplot,wkeep,err,linestyle=1
 oplot,wkeep,-err,linestyle=1
endfor

if(ips ne 0) then device,/close


iprof=0
for ivar=0,nvar-1 do if(varident(ivar,2) eq 0) then iprof=iprof+1

print,'Enter vertical pressure range to plot retrieved profiles over'
yr = fltarr(2)
read,yr

print,'If plotting clouds, plot particles/gram (0) or part/cm3(1)?'
icloud=0
read,icloud

if(ips eq 0)then begin
       window,1,title='Retrievals'
endif else begin
       if(ips eq 1) then begin
        device,filename='retrievals.ps',encapsulated=0
       endif else begin
        device,filename='retrievals.eps',/encapsulated
       endelse
endelse

itypeA=0
print,'Enter nplotx, nploty : '
nplotx=1
nploty=1
read,nplotx,nploty
!p.multi=[0,nplotx,nploty]

readsvp,nsvp,idsvp,svp

istart = 0
for ivar=0,nvar-1 do begin
  itype = varident(ivar,2)
  case itype of
      0: np = npro
      1: np = 2
      2: np = 1
      3: np = 1
      4: np = 3
      6: np = 2
      8: np = 3
      9: np = 3
      888: np = varparam(ivar,0)
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

     rho=xn
     rho(*)=molwt*press/(8.31*temp)

     if(varident(ivar,0) eq 0) then begin
      plot_io,xn,press,yrange=yr,$
      ytitle='Pressure (atm)',xtitle='Temperature (K)',$
      title='Temperature'
      oplot,xn+errn,press,linestyle=1
      oplot,xn-errn,press,linestyle=1
      oplot,xa,press,linestyle=2
      oplot,xa+erra,press,linestyle=3
      oplot,xa-erra,press,linestyle=3
     endif else begin
      if(varident(ivar,0) gt 0) then begin
       tname = string('Gas : ',varident(ivar,0),varident(ivar,1))
       n=0
       j = varident(ivar,0)
       for i = 0,nsvp-1 do if(j eq idsvp(i))then jgas = i
       if(jgas ge 0) then begin
        p1 = exp(svp(0,jgas)+svp(1,jgas)/temp+svp(2,jgas)*temp + svp(3,jgas)*temp^2)
        satvmr = p1/press
        if(j eq 11) then nh3vmr = satvmr
       endif
      endif else begin
       n = abs(varident(ivar,0))
       if(n le ncont) then begin
        tname = string('Cloud type : ',n)
         if(icloud eq 1) then begin
          xn=xn*rho
          xa=xa*rho
          errn=errn*rho
          erra=erra*rho 
         endif    
       endif else begin
        tname='Para-H2 fraction'
       endelse
      endelse
      if(n gt ncont) then begin
        plot_io,xn,press,yrange=yr,$
        ytitle='Pressure (atm)',xrange=[0,1],xstyle=1,$
        title=tname,xtitle='Abundance'
      endif else begin
        if(n eq 0) then begin
         plot_oo,xn,press,yrange=yr,$
         ytitle='Pressure (atm)',xrange=[1e-10,1],$
         title=tname,xtitle='Abundance'
        endif else begin
         b = max([exp(alog(xn)+errn/xn),exp(alog(xa)+erra/xa)])
;         b = 1e5
         plot_io,xn,press,yrange=yr,$
         ytitle='Pressure (atm)',xrange=[0,b],$
         title=tname,xtitle='Abundance'
        endelse
      endelse
      oplot,exp(alog(xn)+errn/xn),press,linestyle=1
      oplot,exp(alog(xn)-errn/xn),press,linestyle=1
      oplot,xa,press,linestyle=2
      oplot,exp(alog(xa)+erra/xa),press,linestyle=3
      oplot,exp(alog(xa)-erra/xa),press,linestyle=3
      if(jgas ge 0) then oplot,satvmr,press,linestyle=4
     endelse
  endif else begin
   if(itype eq 9) then begin

     if(itypeA eq 0) then begin
      print,'Enter gravitational acceleration (for scale height) : '
      g=0.0
      read,g
      scale = 1e-3*8.31*temp/(molwt*g)
     endif

     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))

     h1=xn(2)
     fsh=xn(1)
     xod=xn(0)

     x1=[h1,h1+errn(2),h1-errn(2)]
     p1a = interpol(press,height,x1)
     p1=p1a(0)
     p1p=p1a(1)
     p1m=p1a(2)

     print,'pressure, plus err, minus err',p1,p1p-p1,p1m-p1

     scale1=fsh*scale
     HEa = interpol(scale1,height,x1)
     HE = HEa(0)

     cloud=fltarr(200)
     pX = alog(press)
     pref = pX(0) + (pX(npro-1)-pX(0))*findgen(200)/199.0
     pref = exp(pref)
     href = height(0)+(height(npro-1)-height(0))*findgen(200)/199.0

     ikeep = where(pref lt p1)
     j=ikeep(0)
     cloud(*)=0.
     cloud(j)=xod
     for i=j+1,199 do cloud(i)=xod*exp(-(href(i)-href(j))/HE)
     
     if(itypeA eq 0) then begin
       print,'Enter optical depth range : '
       xr=fltarr(2)
       read,xr

       plot_io,cloud,pref,yrange=yr,xrange=xr

     endif else oplot,cloud,pref

     itypeA=1         

   endif else print,'Profile type not available'

  endelse

  istart = istart+np

endfor

if(ips ne 0) then device,/close


end
