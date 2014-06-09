; *****************************************************************
; Code to compute trapezoidal integation points for an N-point disc
; average calculation
;
; Pat Irwin	4/6/13
;
; *****************************************************************

print,'Enter required number of zenith angles (nstep) : '
nstep=1
read,nstep

mu = 1.0 - findgen(nstep)/float(nstep)
dmu = 1.0/float(nstep)


wt = mu*dmu
wt(0)=0.5*wt(0)

wt=wt*2.
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
tmp(3,*)=tmp(2,*)
tmp(4,*)=180.
tmp(5,*)=wt

print,'Block to paste into spx file (remember to set nav = nstep)'
print,tmp
 
end
