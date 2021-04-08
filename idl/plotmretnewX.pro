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

inpname = strcompress(filename + '.inp',/REMOVE_ALL)
ispace=0
openr,1,inpname
readf,1,ispace
close,1

retname = strcompress(filename + '.mre',/REMOVE_ALL)
apname = strcompress(filename + '.apr',/REMOVE_ALL)
prfname = strcompress(filename + '.prf',/REMOVE_ALL)
refname = strcompress(filename + '.ref',/REMOVE_ALL)
aername = 'aerosol.prf'


readprofilenew,prfname,iformA,nplanetA,xlatA,nproA,ngasA,molwtA,idisoA,heightA,$
        pressA,tempA,vmrA
readrefnew,refname,iformB,nplanetB,xlatB,nproB,ngasB,molwtB,idisoB,heightB,$
        pressB,tempB,vmrB


; ******* Read in .prf file to get pressure grid *******
readprfhead4,prfname,npro,press,height,temp,molwt

if(molwt lt 0.)then begin
 print,'Enter mean molecular weight of atmosphere : '
 read,molwt
 molwt=molwt*1e-3
endif
; ******* Read in aerosol.prf file to get ncont *******
readaerprfhead,aername,ncont

; ******* Read in .apr file to get a priori measurement vector *******
readapriori,apname,npro,nvar,varident,varparam,nx,xa,erra,walb

i31=0

openr,1,retname

 nspec = 1
 readf,1,nspec

 print,nspec

 print,'There are ',nspec,' retrievals in total'
 print,'Enter the number you want to inspect (1-nspec)'
 iplot = 0
 read,iplot
 
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
     nskip=4
;     nskip=44
     for i=1,nskip do begin
      readf,1,head
     endfor
     itype = varident(ivar,2)
     print,'itype = ',itype
     case itype of
     -1: np = npro
      0: np = npro
      1: np = 2
      2: np = 1
      3: np = 1
      4: np = 3
      6: np = 2
      7: np = 2
      8: np = 3
      9: np = 3
      10: np = 4
      11: np = 2
      12: np = 3
      13: np = 3
      14: np = 3
      15: np = 3
      16: np = 4
      17: np = 2
      18: np = 2
      19: np = 4
      20: np = 2
      21: np = 2
      22: np = 5
      25: np = long(varparam(ivar,0))
      29: np = 2*npro
      30: np = long(varparam(ivar,0))
      31: np = long(varparam(ivar,0))
      32: np = 3
      37: np = 1
      38: np = 1
      39: np = 1
      40: np = 2
      41: np = 5
      42: np = 3
      45: np = 3
      102: np = 1
      222: np = 8
      223: np = 9
      224: np = 9
      225: np = 11
      226: np = 8
      333: np = 1
      444: np = 2 + long(varparam(ivar,0))
      555: np = 1
      888: np = long(varparam(ivar,0))
      887: np = long(varparam(ivar,0))
      889: np = 1
      999: np = 1
      777: np = 1
      666: np = 1
     endcase
     if(itype eq 444 and varparam(ivar,1) lt 0) then np = 3
     print,'np = ',np
     xdat1 = fltarr(6,np)
     readf,1,xdat1
     print,'istart = ',istart
     for j=0,np-1 do xdat(*,istart+j)=xdat1(2:5,j)

     istart=istart+np

   endfor

 endfor

close,1

if(ispace eq 0) then begin
 print,'Wavenumber range = ',min(ydat(1,*)),max(ydat(1,*))
endif else print,'Wavelength range = ',min(ydat(1,*)),max(ydat(1,*))

print,'OK (Y/N)?'
ans=''
read,ans
if(strupcase(ans) ne 'Y')then begin
 print,'Enter required plot range : '
 xr = fltarr(2)
 read,xr
endif else begin
 xr = fltarr(2)
 xr(0)=min(ydat(1,*))
 xr(1)=max(ydat(1,*))
endelse

irad=0
irefl=0
print,uname
print,strmid(uname,1,8)
if(strmid(uname,1,8) eq 'Radiance') then begin
 irad=1
 print,'Plot radiance (0) or reflectivity (1) : '
 irefl=0
 read,irefl
endif

if(irefl eq 1) then begin
;  if(ispace eq 1) then begin 
;   openr,3,'/home/oxpln98/plan/irwin/baldrick/raddata/solar.dat'
;   for i=1,3 do readf,3,head
;  endif else begin
;   openr,3,'/home/oxpln98/plan/irwin/baldrick/raddata/solWNum.dat'
;   for i=1,5 do readf,3,head
;  endelse
;  head=' '
;  solar=fltarr(2,78)
;  readf,3,solar
;  close,3

;  if(ispace eq 1) then begin
;  input units are W m-2 nm-1 at 1AU. Convert to uW cm-2 um-1'
;   solar(1,*)=solar(1,*)*1e-4*1e3*1e6
;  endif else begin
;  input units are W m-2 (cm-1)-1 at 1AU. Convert to nW cm-2 (cm-1)-1'
;   solar(1,*)=solar(1,*)*1e-4*1e9   

;   tmp=fltarr(78)
;   tmp(*)=solar(0,*)
;   solar(0,*)=reverse(tmp)
;   tmp(*)=solar(1,*)
;   solar(1,*)=reverse(tmp)
;  endelse
 
 
  if(ispace eq 1) then begin
;   ipfile = '/home/oxpln98/plan/irwin/radtrancode/trunk/raddata/sun_spec_echo.dat'
;   print,'Sun file = ',ipfile
;   openr,1,ipfile
;   head=''
;   for i=1,3 do readf,1,head
;   solar = fltarr(2,600)
;   readf,1,solar
;   close,1

   ipfile = '/network/aopp/oxpln98/plan/irwin/gitradtran/radtrancode/raddata/houghtonsolarwl.dat'
   print,'Sun file = ',ipfile
   print,'OK?'
   ans=''
   read,ans
   ans=strupcase(ans)
   if((ans ne 'Y') and (ans ne 'y')) then begin
    print,'Enter new solfile name : '
    read,ipfile
   endif
   openr,1,ipfile
   head=''
   for i=1,4 do begin
    readf,1,head
    print,head
   endfor
   solx=fltarr(2,9000)
   t1=fltarr(2)
   A=''
   i1=-1
   while ~ eof(1) do begin
    readf,1,A
    print,A
    reads,A,t1
    i1=i1+1
    print,i1,t1
    solx(*,i1)=t1(*)
   endwhile
   nl=i1+1
;   solar = fltarr(2,78)
   solar = fltarr(2,nl)
   for i=0,nl-1 do solar(*,i)=solx(*,i)
;   readf,1,solar
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

if(ips eq 0) then begin
  set_plot,'x'
  device,retain=2
  device,decomposed=0
  window,0,title='Spectrum',xsize=600,ysize=800
endif else begin
 set_plot,'ps'
 if(ips eq 1) then begin
   device,filename='spectrum.ps',encapsulated=0 
 endif else device,filename='spectrum.eps',/encapsulated 
endelse

nconv1=ny/ngeom
!p.multi=[0,1,3]
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
  sun = solint*cos(zenkeep*!pi/180.0)/!pi
  spec=spec/sun
  err=err/sun
  specf=specf/sun
  yname = 'Reflectivity'
endif else begin
 yname=uname
 if(irad eq 1) then begin
  if(ispace eq 0) then begin
   yname = 'Radiance (nW cm!E-2!N sr!E-1!N (cm!E-1!N)!E-1!N)'
  endif else begin
   yname = 'Radiance (!4l!xW cm!E-2!N sr!E-1!N !4l!xm!E-1!N)'
  endelse
 endif
endelse

y2=max([spec+err,specf])

y1=0

print,'y1,y2 = ',y1,y2
print,'Enter new y1,y2'
read,y1,y2

print,'chi/n = ',total(((spec-specf)/err)^2)/float(ny)

xname = 'Wavenumbers (cm!E-1!N)'
if(ispace eq 1) then xname = 'Wavelength (!4l!xm)'

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
  sun = solint*cos(zenkeep*!pi/180.0)/!pi
  spec=spec/sun
  err=err/sun
  specf=specf/sun
 endif


 oplot,wkeep,spec
 oplot,wkeep,spec+err,linestyle=1
 oplot,wkeep,spec-err,linestyle=1
 oplot,wkeep,specf,linestyle=2

endfor

;y1=min([spec-err,specf])
y1=min(specf)


plot_io,wkeep,spec,xtitle=xname,$
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
  sun = solint*cos(zenkeep*!pi/180.0)/!pi
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
 sun = solint*cos(zenkeep*!pi/180.0)/!pi
 spec=spec/sun
 err=err/sun
 specf=specf/sun
endif
yname = 'Difference'

y2=max([specf-spec,err],/NAN)
y1=min([specf-spec,-err],/NAN)
;y2=20*median(abs([specf-spec,err]))
;y1=-y2
plot,wkeep,specf-spec,xtitle=xname,$
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
  sun = solint*cos(zenkeep*!pi/180.0)/!pi
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
for ivar=0,nvar-1 do if(varident(ivar,2) le 0) then iprof=iprof+1

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
  print,'ivar, itype',ivar,itype
  print,'varident = ',varident(ivar,*)

  case itype of
     -1: np = npro
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
      16: np = 4
      17: np = 2
      18: np = 2
      19: np = 4
      20: np = 2
      21: np = 2
      22: np = 5
      30: np = long(varparam(ivar,0))
      31: np = long(varparam(ivar,0))
      32: np = 3
      37: np = 1
      38: np = 1
      39: np = 1
      40: np = 2
      41: np = 5
      42: np = 3
      222: np = 8
      223: np = 9
      224: np = 9
      225: np = 11
      226: np = 8
      333: np = 1
      444: np = 2 + varparam(ivar,0)
      555: np = 1
      888: np = varparam(ivar,0)
      887: np = varparam(ivar,0)
      889: np = 1
      999: np = 1
      777: np = 1
      666: np = 1
  endcase
  jgas=-1

  if(itype le 0) then begin
     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))

     rho=xn

;    Calculate atmospheric density g/cm3.
;    press is already in bars (as read in by readprfhead4.pro)
;    molwt is kg/mol
     rho(*)=molwt*press*100./(8.31*temp)

     if(varident(ivar,0) le 0) then begin
      plot_io,xn,press,yrange=yr,$
      ytitle='Pressure (bar)',xtitle='Temperature (K)',$
      title='Temperature'
      oplot,xn+errn,press,linestyle=1
      oplot,xn-errn,press,linestyle=1
      oplot,xa,press,linestyle=2
;      oplot,xa+10,press,psym=1
      oplot,xa+erra,press,linestyle=3
      oplot,xa-erra,press,linestyle=3

      read,ans

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
        ytitle='Pressure (bar)',xrange=[0,1],xstyle=1,$
        title=tname,xtitle='Abundance'
      endif else begin
        if(n eq 0) then begin
         ikp = where(press le yr(0) and press ge yr(1),nkp)
         a1=min([exp(alog(xn(ikp))-errn(ikp)/xn(ikp)),exp(alog(xa(ikp))-erra(ikp)/xa(ikp))])
         a2=max([exp(alog(xn(ikp))+errn(ikp)/xn(ikp)),exp(alog(xa(ikp))+erra(ikp)/xa(ikp))])
         xr=[a1,a2]
;         xr=[1e-7,0.01]
         plot_oo,xn,press,yrange=yr,xrange=xr,$
         ytitle='Pressure (bar)',$
         title=tname,xtitle='Abundance'
        endif else begin
         ikp = where(press le yr(0) and press ge yr(1),nkp)
         b = max([exp(alog(xn(ikp))+errn(ikp)/xn(ikp)),exp(alog(xa(ikp))+erra(ikp)/xa(ikp))])
;         b = 1e5
         plot_io,xn,press,yrange=yr,$
         ytitle='Pressure (bar)',xrange=[0,b],$
         title=tname,xtitle='Abundance'
        endelse
      endelse
      oplot,exp(alog(xn)+errn/xn),press,linestyle=1
      oplot,exp(alog(xn)-errn/xn),press,linestyle=1
      oplot,xa,press,linestyle=2
      oplot,exp(alog(xa)+erra/xa),press,linestyle=3
      oplot,exp(alog(xa)-erra/xa),press,linestyle=3
      if(jgas ge 0) then oplot,satvmr,press,linestyle=4

      read,ans

     endelse
  endif else begin
   if(itype eq 9) then begin

     if(itypeA eq 0) then begin
      print,'Enter gravitational acceleration (for scale height) : '
      g=0.0
      read,g
;     g is m s-2
;     temp is K
;     molwt is kg/mol
;     R is J mol-1 K-1
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
       itypeA=1

     endif else begin

      oplot,cloud,pref,linestyle=ivar
      print,'Overplotting'
      print,'max, min = ',max(cloud),min(cloud)
      read,ans
     endelse

;     read,ans

   endif else begin

    if(itype ge 12 and itype le 15) then begin

     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))

;     xod=exp(xn(0))
;     p1=exp(xn(1))
;     wid=exp(xn(2))

     xod=xn(0)
     wid=xn(2)
     if (itype le 13) then p1=xn(1) else h1=xn(1)
     if (itype le 13) then lp1 = alog(p1)

     cloud=fltarr(200)
     pX = alog(press)
     hX = height
     tX = temp

     pref = pX(0) + (pX(npro-1)-pX(0))*findgen(200)/199.0
     href = hX(0) + (hX(npro-1)-hX(0))*findgen(200)/199.0
     tref = interpol(temp,height,href)

;    pressure already in bars
;    molwt in kg/mol
;    Calculate density in g/cm3
     rho=molwt*exp(pref)*100/(8.31*tref)


     case itype of 
      12:  cloud = xod*exp(-((pref-lp1)/wid)^2)
      13:  cloud = (xod*wid^2)/(wid^2 + (pref-lp1)^2)
      14:  cloud = (xod/wid)*exp(-((href-h1)/wid)^2)
      15:  cloud = (xod*wid)/(wid^2 + (href-h1)^2)
     endcase


     for i=0,npro-1 do print,i,height(i),cloud(i),cloud(i)*rho(i)

     if(icloud eq 1) then begin
      cloud=cloud*rho
     endif

     print,'Max abundance of this cloud = ',max(cloud)

     if(itypeA eq 0) then begin
       print,'Enter max abundance to plot : '
       x1=1.
       read,x1
       xr=[0,x1]

       plot_io,cloud,exp(pref),yrange=yr,xrange=xr

     endif else oplot,cloud,exp(pref)

;     itypeA=1         

    endif else begin
     if(itype eq 3) then begin
      jgas=-1
      for igas=0,ngasA-1 do begin
       if(varident(ivar,0) eq idisoA(0,igas) and varident(ivar,1) eq idisoA(1,igas))then begin
        jgas=igas
       endif
      endfor

      if(jgas eq -1 and varident(ivar,0) eq 0)then begin
       plot_io,tempA,pressA,yrange=yr
       oplot,tempB,pressB,linestyle=1
      endif else begin
       plot_oo,vmrA(*,jgas),pressA,yrange=yr
       oplot,vmrB(*,jgas),pressB,linestyle=1
      endelse

      read,ans

     endif else begin 
      if(itype ne 16 and itype ne 444 and itype ne 222) then begin
        print,'Profile type not available'
      endif
     endelse
    endelse
   endelse

  endelse

  if(itype eq 30) then begin
   nlevel = long(varparam(ivar,1))  
   nlong = np/nlevel
   long = 360.*findgen(nlong+1)/float(nlong)

   xn = fltarr(np)
   xa = fltarr(np)
   errn = xn
   erra = xa
   xa(*)=xdat(0,istart:(istart+np-1))
   erra(*)=xdat(1,istart:(istart+np-1))
   xn(*)=xdat(2,istart:(istart+np-1))
   errn(*)=xdat(3,istart:(istart+np-1))
 
   p1 = fltarr(nlevel)
   p1(*)=varparam(ivar,2:2+nlevel-1)

   xall = fltarr(nlong+1,nlevel)
   ximprove = xall
   for ilong = 0,nlong-1 do for ilevel = 0,nlevel-1 do begin
    k = ilong*nlevel+ilevel
    xall(ilong,ilevel)=xn(k)
    ximprove(ilong,ilevel)=1. - errn(k)/erra(k)
   endfor
   xall(nlong,*)=xall(0,*)
   window,3
   !p.multi=[0,1,2]

   contour,xall,long,p1,xtitle='Longitude',ytitle='pressure',$
    xstyle=1,/follow,/ylog,yrange=[max(p1),min(p1)],nlevel=10,title='Xn'

   contour,ximprove,long,p1,xtitle='Longitude',ytitle='pressure',$
    xstyle=1,/follow,/ylog,yrange=[max(p1),min(p1)],nlevel=10,$
    title='Improvement factor'

   print,'Enter pressure level to  contour plot : '
   px = 1.0
   read,px

   ikeep = where(px gt p1)
   ilev=ikeep(0)

   print,'ilev,plev = ',ilev,p1(ilev)

   nlat=7
   lat = -90. + 30.*findgen(nlat)
   tmap = fltarr(nlong+1,7)
   xlong = fltarr(nlong)
   for ilong = 0,nlong-1 do xlong(ilong) = xn(ilong*nlevel+ilev)
   sum=total(xlong)/float(nlong)
   for ilong=0,nlong-1 do begin
    for ilat=0,nlat-1 do begin
     tmap(ilong,ilat)=sum+(xlong(ilong)-sum)*cos(lat(ilat)*!pi/180.)
    endfor
   endfor
   tmap(nlong,*)=tmap(0,*)

   window,4
   !p.multi=0
   contour,tmap,long,lat,yrange=[-90,90.],xrange=[0,360.],$
    xstyle=1,ystyle=1,nlevel=10,/follow
   wset,1
   !p.multi=[0,nplotx,nploty]
  endif

  if(itype eq 31) then begin
   nlong = long(varparam(ivar,0))  
   long = 360.*findgen(nlong+1)/float(nlong)

   xn = fltarr(np)
   xa = fltarr(np)
   errn = xn
   erra = xa
   xa(*)=xdat(0,istart:(istart+np-1))
   erra(*)=xdat(1,istart:(istart+np-1))
   xn(*)=xdat(2,istart:(istart+np-1))
   errn(*)=xdat(3,istart:(istart+np-1))
 
   xall = fltarr(nlong+1)
   exalln = xall
   exalla = xall
  
   xall(0:nlong-1)=xn(*)
   xall(nlong)=xn(0)
   exalln(0:nlong-1)=errn(*)
   exalln(nlong)=errn(0)
   exalla(0:nlong-1)=erra(*)
   exalla(nlong)=erra(0)

   f = fltarr(nlong+1)
   fa = f
   fa(*)=1.0

   vref = [1e-5,9.9e-8,1e-5,1e-5]
   idref = [1,2,5,6]

   if(i31 eq 0) then begin
     window,3
     !p.multi=0
     plot_io,long,xall*vref(i31),xrange=[0,360],xstyle=1,yrange=[1e-10,1],$
      ystyle=1
   endif

   oplot,long,xall*vref(i31),linestyle=0

   f(*)=exalln/xall

   f1 = sqrt(1.0/(1.0/f^2 - 1/fa^2))
   oplot,long,xall*vref(i31)*(1.0+f1),linestyle=1
   oplot,long,xall*vref(i31)/(1.0+f1),linestyle=1
   print,'Gas ID = ',idref(i31)

   ans=''
   read,ans

   i31=i31+1

  endif


  if(itype eq 444) then begin

     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))

     plot,walb,xn(2:np-1),xtitle='Wavelength',ytitle='Imaginary RI'
     oplot,walb,xn(2:np-1)+errn(2:np-1),linestyle=1
     oplot,walb,xn(2:np-1)-errn(2:np-1),linestyle=1
     oplot,walb,xa(2:np-1),linestyle=2
     oplot,walb,xa(2:np-1)-erra(2:np-1),linestyle=2
     oplot,walb,xa(2:np-1)+erra(2:np-1),linestyle=2

     print,'Apriori radius and error : ',xa(0),erra(0)
     print,'Apriori variance and error : ',xa(1),erra(1)
     print,'Fitted radius and error : ',xn(0),errn(0)
     print,'Fitted variance and error : ',xn(1),errn(1)

  endif


  if(itype eq 222) then begin

     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))

     plot,xn(2:np-1),ytitle='Sromovsky Cloud Parameters'
     oplot,xn(2:np-1)+errn(2:np-1),linestyle=1
     oplot,xn(2:np-1)-errn(2:np-1),linestyle=1
     oplot,xa(2:np-1),linestyle=2
     oplot,xa(2:np-1)-erra(2:np-1),linestyle=2
     oplot,xa(2:np-1)+erra(2:np-1),linestyle=2

  endif
 
  

  if(itype eq 16) then begin
     xn = fltarr(np)
     xa = fltarr(np)
     errn = xn
     erra = xa
     xa(*)=xdat(0,istart:(istart+np-1))
     erra(*)=xdat(1,istart:(istart+np-1))
     xn(*)=xdat(2,istart:(istart+np-1))
     errn(*)=xdat(3,istart:(istart+np-1))

     xxa=fltarr(npro)
     xxa(*)=xa(0)
     
     pknee = xa(1)
     hknee = interpol(height,press,pknee)
     ikeep = where(press lt pknee)
     jref = ikeep(0)-1

     for i=jref,0,-1 do begin
      delh = hknee-height(i)
      xxa(i)=xa(0)+xa(2)*delh
     endfor
     for i=jref+1,npro-1 do begin
      delh = height(i)-hknee
      xxa(i)=xa(0)+xa(3)*delh
     endfor

     plot_io,xxa,press,yrange=yr

     pknee = xn(1)
     hknee = interpol(height,press,pknee)
     ikeep = where(press lt pknee)
     jref = ikeep(0)-1

     for i=jref,0,-1 do begin
      delh = hknee-height(i)
      xxa(i)=xn(0)+xn(2)*delh
     endfor
     for i=jref+1,npro-1 do begin
      delh = height(i)-hknee
      xxa(i)=xn(0)+xn(3)*delh
     endfor

     oplot,xxa,press,linestyle=1

  endif

  istart = istart+np

endfor

if(ips ne 0) then device,/close


end
