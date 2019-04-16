; *************************************************************************************************
; Procedure to read in internal flux output file, ish.dat, and construct internal radiation field.
;
; Pat Irwin     8/10/17
;
; *************************************************************************************************

openr,1,'ish.dat'
nlayer=1
nmu=0
nwave=1
nf=1
sol_ang=0.

readf,1,nlayer
readf,1,nmu
readf,1,nwave
readf,1,nf
readf,1,sol_ang

mu0 = cos(sol_ang*!pi/180.)
wave=fltarr(nwave)
solar=wave
radground=wave
galb=wave
basep=fltarr(nlayer)
baseh=basep

readf,1,wave
readf,1,basep
readf,1,baseh

basep=reverse(basep)
baseh=reverse(baseh)

height = [2*baseh(0),baseh]
press = [0.,basep] 
; xgplf(nwave,nmu,nlayer,nf+1) Internal radiances (downwards) at bottom of each layer
; xgmif(nwave,nmu,nlayer,nf+1) Internal radiances (upwards) at top of each layer
; xgplfb(nwave,nmu,nlayer,nf+1) Internal radiances (downwards) at bottom of each layer for non-scattering conditions (to get direct solar)
; xgmifb(nwave,nmu,nlayer,nf+1) Internal radiances (upwards) at top of each layer for non-scattering conditions (direct solar)
; xgplfa(nwave,nmu,nlayer) Internal radiances (downwards) at bottom of each layer for scattering thermal emission only
; xgmifa(nwave,nmu,nlayer) Internal radiances (upwards) at top of each layer for scattering thermal emission only
; xgplfc(nwave,nmu,nlayer) Internal radiances (downwards) at bottom of each layer for non-scattering thermal emission only
; xgmifc(nwave,nmu,nlayer) Internal radiances (upwards) at top of each layer for non-scattering thermal emission only

xgplf=fltarr(nwave,nmu,nlayer,nf+1)
xgmif=fltarr(nwave,nmu,nlayer,nf+1)
xgplfb=fltarr(nwave,nmu,nlayer,nf+1)
xgmifb=fltarr(nwave,nmu,nlayer,nf+1)
xgplfa=fltarr(nwave,nmu,nlayer)
xgmifa=fltarr(nwave,nmu,nlayer)
xgplfc=fltarr(nwave,nmu,nlayer)
xgmifc=fltarr(nwave,nmu,nlayer)

directsol = fltarr(nwave,2*nmu,nlayer+1,360)
diffusescat = fltarr(nwave,2*nmu,nlayer+1,360)
thermal = fltarr(nwave,2*nmu,nlayer+1)
thermalb = fltarr(nwave,2*nmu,nlayer+1)
phi = findgen(360)*!pi/180.
drad = fltarr(360)

directsol(*)=0.
diffusescat(*)=0.
thermal(*)=0.
thermalb(*)=0.

;fluxes
fplf = fltarr(nwave,nlayer+1)
fplfb = fltarr(nwave,nlayer+1)
fmif = fltarr(nwave,nlayer+1)
fmifb = fltarr(nwave,nlayer+1)

tmp=fltarr(1,nf+1)

; Now read in radiance arrays
x=1.
for iwave=0,nwave-1 do begin
 readf,1,x
 solar(iwave)=x
 readf,1,x
 radground(iwave)=x
 readf,1,x
 galb(iwave)=x
 for imu=0,nmu-1 do for ilay=0,nlayer-1 do begin
  readf,1,tmp
  xgplf(iwave,imu,ilay,*)=tmp(0,*)
  readf,1,tmp
  xgmif(iwave,imu,ilay,*)=tmp(0,*)
  readf,1,tmp
  xgplfb(iwave,imu,ilay,*)=tmp(0,*)
  readf,1,tmp
  xgmifb(iwave,imu,ilay,*)=tmp(0,*)
  readf,1,x
  xgplfa(iwave,imu,ilay)=x
  readf,1,x
  xgmifa(iwave,imu,ilay)=x
  readf,1,x
  xgplfc(iwave,imu,ilay)=x
  readf,1,x
  xgmifc(iwave,imu,ilay)=x
 endfor
endfor

close,1

print,'Enter calculation run name : '
ipfile=''
read,ipfile
ipfile=ipfile+'.sca'
openr,1,ipfile
head=''
readf,1,head
readf,1,head
xt = fltarr(2,nmu)
readf,1,xt
close,1


mu=fltarr(nmu)
wt=fltarr(nmu)
mu(*)=xt(0,*)
wt(*)=xt(1,*)

iord=findgen(nmu)

ix = interpol(iord,mu,mu0)
imu0 = long(ix+0.001)

dimu = ix-imu0
print,'sol_ang, mu0 = ',sol_ang,mu0
print,'imu0,mu(imu0),dimu',imu0,mu(imu0),dimu



;Calculate fluxes
xnorm=2*!pi
for iwave=0,nwave-1 do for ilayer=0,nlayer do begin
 case ilayer of 
 0: begin
         fplf(iwave,ilayer) = xnorm*mu(imu0)*wt(imu0)*solar(iwave)/(2*!pi*wt(imu0))
         fplfb(iwave,ilayer) = xnorm*mu(imu0)*wt(imu0)*solar(iwave)/(2*!pi*wt(imu0))
         fmif(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgmif(iwave,*,ilayer,0))
         fmifb(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgmifb(iwave,*,ilayer,0))
    end

 nlayer: begin
          fplf(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgplf(iwave,*,ilayer-1,0))
          fplfb(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgplfb(iwave,*,ilayer-1,0))
          fmif(iwave,ilayer) = fplf(iwave,ilayer)*galb(iwave)+radground(iwave)
          fmifb(iwave,ilayer) = fplfb(iwave,ilayer)*galb(iwave)+radground(iwave)
         end

 else: begin
         fplf(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgplf(iwave,*,ilayer-1,0))
         fplfb(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgplfb(iwave,*,ilayer-1,0))
         fmif(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgmif(iwave,*,ilayer,0))
         fmifb(iwave,ilayer) = xnorm*total(mu(*)*wt(*)*xgmifb(iwave,*,ilayer,0))
       end
 endcase
endfor


; Need to regrid and reorder radiances.
; Combine radiances so that the upwards and downwards relate to same pressure levels at top and bottom of
; each layer
; Also reorder so that directions of 2*nmu cos(theta) ordinates are: up, slightly up, side, slightly down, down

; Do thermal first
for ilayer=0,nlayer do begin
 case ilayer of
  0:begin
        for iwave=0,nwave-1 do begin
         for imu=0,nmu-1 do begin
          thermal(iwave,nmu-1-imu,ilayer) = xgmifa(iwave,imu,ilayer) 
          thermal(iwave,nmu+imu,ilayer) = 0. 
          thermalb(iwave,nmu-1-imu,ilayer) = xgmifc(iwave,imu,ilayer) 
          thermal(iwave,nmu+imu,ilayer) = 0. 
          
          for ic=0,nf do begin
           drad(*) = xgmif(iwave,imu,ilayer,ic)*cos(ic*phi)
           if(ic eq 0) then drad(*)=drad(*)-xgmifa(iwave,imu,ilayer) else drad(*)=drad(*)*2.
           diffusescat(iwave,nmu-1-imu,ilayer,*)=diffusescat(iwave,nmu-1-imu,ilayer,*)+drad(*)
          endfor

         endfor
     
         for ic=0,nf do begin
           drad(*) = solar(iwave)/(2*!pi*wt(imu0))*cos(ic*phi)
           if(ic gt 0) then drad(*)=drad(*)*2.
;           print,ilayer,ic,drad(*)
           diffusescat(iwave,nmu+imu0,ilayer,*)=diffusescat(iwave,nmu+imu0,ilayer,*)+drad(*)
         endfor
         directsol(iwave,nmu+imu0,ilayer,*) = diffusescat(iwave,nmu+imu0,ilayer,*)

        endfor
    end

  nlayer:begin
	   for iwave=0,nwave-1 do begin
            for imu=0,nmu-1 do begin
             thermal(iwave,nmu-1-imu,ilayer) = radground(iwave) 
             thermal(iwave,nmu+imu,ilayer) = xgplfa(iwave,imu,ilayer-1) 
             thermalb(iwave,nmu-1-imu,ilayer) = radground(iwave) 
             thermalb(iwave,nmu+imu,ilayer) = xgplfc(iwave,imu,ilayer-1) 
             diffusescat(iwave,nmu-1-imu,ilayer,*) = fplf(iwave,ilayer)*galb(iwave)/!pi
; don't need?           directsol(iwave,nmu-1-imu,ilayer,*) = fplfb(iwave,ilayer)*galb(iwave)/!pi
 
             for ic=0,nf do begin

              drad(*) = xgplf(iwave,imu,ilayer-1,ic)*cos(ic*phi)
              if(ic eq 0) then drad(*)=drad(*)-xgplfa(iwave,imu,ilayer-1) else drad(*)=drad(*)*2.
              diffusescat(iwave,nmu+imu,ilayer,*)=diffusescat(iwave,nmu+imu,ilayer,*)+drad(*)

              if(imu eq imu0) then begin
               drad(*) = xgplfb(iwave,imu,ilayer-1,ic)*cos(ic*phi)
               if(ic eq 0) then drad(*)=drad(*)-xgplfc(iwave,imu,ilayer-1) else drad(*)=drad(*)*2.
               directsol(iwave,nmu+imu,ilayer,*)=directsol(iwave,nmu+imu,ilayer,*)+drad(*)
              endif

             endfor
            endfor
           endfor
         end

  else:begin
        for iwave=0,nwave-1 do begin
         for imu=0,nmu-1 do begin

          thermal(iwave,nmu-1-imu,ilayer) = xgmifa(iwave,imu,ilayer) 
          thermal(iwave,nmu+imu,ilayer) = xgplfa(iwave,imu,ilayer-1) 
          thermalb(iwave,nmu-1-imu,ilayer) = xgmifc(iwave,imu,ilayer) 
          thermalb(iwave,nmu+imu,ilayer) = xgplfc(iwave,imu,ilayer-1) 

          for ic=0,nf do begin

           drad(*) = xgmif(iwave,imu,ilayer,ic)*cos(ic*phi)
           if(ic eq 0) then drad(*) = drad(*)-xgmifa(iwave,imu,ilayer) else drad(*)=drad(*)*2.
           diffusescat(iwave,nmu-1-imu,ilayer,*) = diffusescat(iwave,nmu-1-imu,ilayer,*) + drad(*)

           drad(*) = xgplf(iwave,imu,ilayer-1,ic)*cos(ic*phi)
           if(ic eq 0) then drad(*)=drad(*)-xgplfa(iwave,imu,ilayer-1) else drad(*)=drad(*)*2.
           diffusescat(iwave,nmu+imu,ilayer,*) = diffusescat(iwave,nmu+imu,ilayer,*) + drad(*)

;           if(imu eq imu0 and ilayer eq 1) then print,ilayer,ic,drad(*)

; don't     drad(*) = xgmifb(iwave,imu,ilayer,ic)*cos(ic*phi)
; need      if(ic eq 0) then drad(*)=drad(*)-xgmifc(iwave,imu,ilayer) else drad(*)=drad(*)*2.
;           directsol(iwave,nmu-1-imu,ilayer,*) = directsol(iwave,nmu-1-imu,ilayer,*) + drad(*)

           if(imu eq imu0) then begin
            drad(*) = xgplfb(iwave,imu,ilayer-1,ic)*cos(ic*phi)
            if(ic eq 0) then drad(*)=drad(*)-xgplfc(iwave,imu,ilayer-1) else drad(*)=drad(*)*2.
            directsol(iwave,nmu+imu,ilayer,*) = directsol(iwave,nmu+imu,ilayer,*) + drad(*)
           endif

          endfor
         endfor
        endfor
       end
 endcase

endfor

; Assume top layer has effectively zero opacity and set flux down at top of top layer to be flux down at bottom of top layer
;for iwave=0,nwave-1 do diffusescat(iwave,nmu+imu0,0,*)=diffusescat(iwave,nmu+imu0,1,*)


print,'Enter level to plot scattered radiation in range 0 to ',nlayer

read,ilayer
iwave=0
plotdiffuse1,diffusescat,directsol,mu,nwave,nmu,nlayer,iwave,ilayer


end
