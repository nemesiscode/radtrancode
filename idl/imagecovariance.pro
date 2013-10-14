; ********************************************************************
; Program to graphically display the diagnostic information written to
; the .cov file by Nemesis for the case of a single retrieval
;
; Pat Irwin	10/9/04
;
; ********************************************************************

print,'Enter input filename ' 
filename=''
read,filename

openr,1,filename
npro=1
nvar=1
readf,1,npro,nvar
itmp = lonarr(3)
tmp = fltarr(5)
varident = lonarr(nvar,3)
varparam = fltarr(nvar,5)

for i=0,nvar-1 do begin
 readf,1,itmp
 readf,1,tmp
 varident(i,*)=itmp(*)
 varparam(i,*)=tmp(*)
endfor

nx=1
ny=1
readf,1,nx,ny

sa = fltarr(nx,nx)
st = sa
sm = sa
sn = sa

tmp = fltarr(nx)
for i=0,nx-1 do begin
 readf,1,tmp
 sa(i,*)=tmp(*)
 readf,1,tmp
 sm(i,*)=tmp(*)
 readf,1,tmp
 sn(i,*)=tmp(*)
 readf,1,tmp
 st(i,*)=tmp(*)
endfor

aa = fltarr(nx,nx)
readf,1,aa

dd = fltarr(nx,ny)
readf,1,dd

kk = fltarr(ny,nx)
readf,1,kk

kt = transpose(kk)

se1 = fltarr(ny)
readf,1,se1

close,1

nplotx = 2
nploty = 2
ict = 33
ips = 0

x1=0.0
x2=0.0
imagematrix,sa,x1,x2,'sa',1,nplotx,1,nploty,ips,ict
imagematrix,sm,x1,x2,'sm',2,nplotx,1,nploty,ips,ict
imagematrix,sn,x1,x2,'sn',1,nplotx,2,nploty,ips,ict
imagematrix,st,x1,x2,'st',2,nplotx,2,nploty,ips,ict

ans=''
read,ans

!p.multi=0
xn = findgen(nx)
ea = fltarr(nx)
em = ea
en = ea
et = ea
for i=0,nx-1 do begin
 ea(i)=sqrt(sa(i,i))
 en(i)=sqrt(sn(i,i))
 em(i)=sqrt(sm(i,i))
 et(i)=sqrt(st(i,i))
endfor

x2 = max([ea,en,em,et])
!p.position=0
!p.multi=0
plot,et,xn,xtitle='Error',ytitle='xn',xrange=[0,x2]
oplot,ea,xn,linestyle=1
oplot,en,xn,linestyle=2
oplot,em,xn,linestyle=3
print,'linestyle 0 = e_total'
print,'linestyle 1 = e_apriori'
print,'linestyle 2 = e_smooth'
print,'linestyle 3 = e_meas'
read,ans

ca = sa
cm = sm
cn = sn
ct = st

for i=0,nx-1 do begin
 for j=0,nx-1 do begin
  ca(i,j)=sa(i,j)/sqrt(sa(i,i)*sa(j,j))
  cm(i,j)=sm(i,j)/sqrt(sm(i,i)*sm(j,j))
  cn(i,j)=sn(i,j)/sqrt(sn(i,i)*sn(j,j))
  ct(i,j)=st(i,j)/sqrt(st(i,i)*st(j,j))
 endfor
endfor

x1 = -1.
x2 = 1.0
imagematrix,ca,x1,x2,'ca',1,nplotx,1,nploty,ips,ict
imagematrix,cm,x1,x2,'cm',2,nplotx,1,nploty,ips,ict
imagematrix,cn,x1,x2,'cn',1,nplotx,2,nploty,ips,ict
imagematrix,ct,x1,x2,'ct',2,nplotx,2,nploty,ips,ict

openw,22,'correlation.txt'
printf,22,ct
close,22

read,ans

nplotx=1
nploty=1
x1 = 0.0
x2 = 0.0
imagematrix,aa,x1,x2,'aa',1,nplotx,1,nploty,ips,ict

read,ans
!p.position=0

!p.multi=[0,3,1]
x1 = min(aa)
x2 = max(aa)
xn = findgen(nx)
plot,aa(*,0),xn,xrange=[x1,x2]
for i=1,nx-1 do oplot,aa(*,i),xn

stot = fltarr(nx)
for i=0,nx-1 do stot(i)=total(aa(*,i))

fwhm = fltarr(nx)
for i=0,nx-1 do begin
j=i
amax = aa(i,i)
st1: j=j+1
     if(j eq nx)then goto,st2
     at = aa(j,i)
     if(at gt 0.5*amax)then goto,st1

st2: k=i
     
st3: k=k-1
     if(k lt 0) then goto,st4
     at = aa(k,i)
     if(at gt 0.5*amax)then goto,st3

st4: del = j-k
fwhm(i)=del
endfor

plot,stot,xn,title='Averaging kernel area (i.e. total)'

plot,fwhm,xn,title='FWHM'

read,ans
!p.multi=0

x1=0
x2=0 
imagematrix,dd,x1,x2,'dd',1,nplotx,1,nploty,ips,ict

read,ans
!p.multi=0
x1 = min(dd)
x2 = max(dd)

plot,dd(*,0),xn,xrange=[x1,x2],title='Contribution functions'
for i=1,ny-1 do oplot,dd(*,i),xn

read,ans
!p.multi=0

x1=0
x2=0 
imagematrix,kk,x1,x2,'kk',1,nplotx,1,nploty,ips,ict

read,ans

se = fltarr(ny,ny)
se(*)=0.0
for i=0,ny-1 do se(i,i)=se1(i)

nplotx=2
nploty=2

s1 = sa#kt
s1 = kk#s1

x1=0
x2=0
imagematrix,s1,x1,x2,'kk#sa#kt',1,nplotx,1,nploty,ips,ict
imagematrix,se,x1,x2,'se',2,nplotx,1,nploty,ips,ict
imagematrix,s1+se,x1,x2,'kk#sa#kt+ se',1,nplotx,2,nploty,ips,ict

read,ans
!p.position=0
!p.multi=0

se2 = fltarr(ny)
for i=0,ny-1 do se2(i)=s1(i,i)

x1 = min([se1,se2])
x2 = max([se1,se2])

dum=indgen(ny)
plot_io,dum,se1,yrange=[x1,x2],title='Solid = measurement, dots: kk#sa#kt'
oplot,dum,se2,linestyle=1

end
