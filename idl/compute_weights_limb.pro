; *****************************************************************
; Code to compute trapezoidal integation points for an N-point disc
; average calculation, paying particular attention to the limb region 
; assuming that IPZEN=2 and thus that the emission angles are defined at
; the TOP of the atmosphere.
;
; User has to specify how many points with paths intersecting the solid
; surface of the planet and how many limb paths.
;
; Code needs the .drv file computed for a nadir path in order to get the
; base heights of the layers.
;
; Pat Irwin	22/5/15
;
; *****************************************************************

print,'Enter radius of planet : '
radius=1
read,radius

print,'Enter name of a .drv file for a nadir path (to get layer base heights)'
drvname=''
read,drvname

subreaddrv4,drvname,nlay,npath,ngas,ncont,idgas,isogas,iprogas,$
baseh,delh,basep,baset,totam,press,temp,dopp,ppgas,amgas,cloud


print,'Enter required number of zenith angles from centre to edge of surface : '
nstep1=1
read,nstep1

print,'Enter number of steps from edge of surface to space : '
nstep2=1
read,nstep2


; Find edge of planet's surface
Rtop = baseh(nlay-1)+delh(nlay-1)
thetaedge = asin(radius/(radius+Rtop))
muedge=cos(thetaedge)

;find mu for all base of all layers. 
thetalay=asin((radius+baseh)/(radius+Rtop))
mulay=cos(thetalay)

; find mu for heights just above the base of each layer. These are to 
; guard agains numerical inconsistencies between this code and Nemesis to
; ensure that Nemesis uses the same layers as is intended here
mulayX=mulay
for i=0,nlay-1 do begin
 rr = baseh(i)+0.1*delh(i)
 if(radius+rr le radius+Rtop)then begin
  theta = asin((radius+rr)/(radius+Rtop))
  mulayX(i)=cos(theta)
 endif else mulayX(i)=0.
endfor

; mu values for paths intersecting surface
mu1 = 1.0 - (1.0-muedge)*findgen(nstep1)/float(nstep1)

; mu values for limb paths
mu2=fltarr(nstep2)
mu2X=fltarr(nstep2)
for i=0,nstep2-1 do begin
 xmu = muedge - muedge*i/float(nstep2-1)
 ikeep = where(mulay ge xmu,nmu)
 jlay=ikeep(nmu-1)
 mu2(i)=mulay(jlay) 
 mu2X(i)=mulayX(jlay) 
endfor

; Strip out any repeated values of mu2
mu3=fltarr(nstep2)
mu3(0)=mu2(0)
nn=0
for i=1,nstep2-1 do begin
 if(mu2(i) ne mu3(nn))then begin
  nn=nn+1
  mu3(nn)=mu2(i)
 endif
endfor

nstep2=nn+1
mu2=fltarr(nstep2)
mu2(*)=mu3(0:nn)


; final set of mu values
mu=[mu1,mu2]
nstep=nstep1+nstep2

;Now compute weights for disc-integration using Trapezium rule
wt=fltarr(nstep)
for i=0,nstep-1 do begin
 if(i eq 0) then begin
  delmu1=mu(i)-mu(i+1)
  wt(i)=mu(i)*delmu1
 endif else begin
  delmu1=mu(i-1)-mu(i)
  if(i lt nstep-1) then begin
   delmu2=mu(i)-mu(i+1)
  endif else delmu2=mu(i)
  wt(i)=mu(i)*(delmu1+delmu2)
 endelse  
endfor

print,'  i     mu(i)    wt(i)'
for i=0,nstep-1 do print,i,mu(i),wt(i)
print,'Sum of weights = ',total(wt)

; ********************************************************
; Now part to output FOV averaging block to paste into an spx file:

print,'Enter latitude and longitude of point : '
xlat=1.
xlon=1.
read,xlat,xlon

; set viewing and solar zenith angles to be the same as acos(mu).

tmp=fltarr(6,nstep)
tmp(0,*)=xlat
tmp(1,*)=xlon
tmp(2,*)=acos(mu)*180./!pi
; For tangent paths make sure the zenith angles are very slightly bigger to 
; ensure we that the bottom layer assigned by Nemesis is as intended. 
for i=nstep1,nstep-1 do tmp(2,i)=acos(mu2X(i-nstep1))*180./!pi
tmp(3,*)=tmp(2,*)
tmp(4,*)=180.
tmp(5,*)=wt

print,'Block to paste into spx file (remember to set nav = nstep)'
print,'nstep = ',nstep
print,tmp
 
end
