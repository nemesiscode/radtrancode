npoint=121
k = findgen(npoint)/100.

Imeanline=k
Imeandisc=2/(1+2*k)

for i=0,npoint-1 do Imeanline(i)=calc_int(k(i))

xx = [0.5,0.5]
yy = [0,2]
Iref=k
Iref(*)=1.0
!p.multi=[0,1,2]
plot,k,Imeandisc,xtitle='minnaert-k',ytitle='I_mean/I0',$
title='Imean/I0 as function of k'
oplot,k,Imeanline,linestyle=1
oplot,xx,yy,linestyle=2
oplot,k,Iref,linestyle=2

thet=45.0
mu1 = cos(thet*!pi/180.0)

xbeta = mu1^(1.0-2*k)

plot,k,Imeandisc*xbeta,xtitle='minnaert-k',ytitle='I_mean/I(45)',yrange=[0,2],$
title='Imean/I(45) as function of k'
oplot,k,Imeanline*xbeta,linestyle=1
oplot,xx,yy,linestyle=2
oplot,k,Iref,linestyle=2

mu = [1.00000000000000,0.919533908166459,0.738773865105505,0.477924949810444,0.165278957666387]
wt = [2.222222222222220E-002,0.133305990851069,0.224889342063117,0.292042683679684,0.327539761183898]

thet = acos(mu)
thet = [thet,0]

dwtdisc = 2*mu*wt

delthet = fltarr(5)
for i=0,3 do delthet(i)=thet(i+1)-thet(i)
delthet(4)=0.5*!pi-thet(4)

dwtline = dwtdisc
dwtline(0)=0.5*mu(0)*delthet(0)
for i=1,4 do dwtline(i)=0.5*mu(i)*(delthet(i-1)+delthet(i))


; numerical re-normalisation

dwtline = dwtline/total(dwtline)
dwtdisc = dwtdisc/total(dwtdisc)

print,'Quadrature zenith angles (degree) : '
print,thet(0:4)*180/!pi
print,'Disc-integration weights : '
print,dwtdisc
print,'Line-integration weights : '
print,dwtline

window,1
plot,thet(0:4)*180/!pi,dwtdisc,xtitle='zenith angle',ytitle='Weight',$
title='Weights for disc and line-averaging'
oplot,thet(0:4)*180/!pi,dwtline,linestyle=1

plot,mu,dwtdisc,xtitle='Cos(zenith)',ytitle='Weight',$
title='Weights for disc and line-averaging'

oplot,mu,dwtline,linestyle=1

openw,1,'lineminnaert.txt'
printf,1,npoint
tmp=fltarr(2,npoint)
tmp(0,*)=k(*)
tmp(1,*)=Imeanline(*)
printf,1,tmp
close,1

end
