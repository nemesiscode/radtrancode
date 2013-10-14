pro imagematrix,kkin,kmin,kmax,tname,iplotx,nplotx,iploty,nploty,ips,ict

xwid = !d.x_vsize
ywid = !d.y_vsize

; print,iplotx,iploty

x0 = float(iplotx-1)/float(nplotx)
y0 = float(nploty - iploty)/float(nploty)
dx = 1.0/float(nplotx)
dy = 1.0/float(nploty)

x1 = x0+dx
y1 = y0+dy
; print,'plotbars : ',x0,y0,x1,y1

if(iplotx eq 1) then begin
  lm = 0.3*dx
  rm = 0.1*dx
endif else begin
  rm = 0.3*dx
  lm = 0.1*dx
endelse

bm = 0.3*dy
tm = 0.2*dy

itmp = size(kkin)
nx = itmp(1)
ny = itmp(2)
x=[1,nx]
y=[1,ny]

!p.position = [x0+lm,y0+bm,x1-rm,y1-tm]
if(iplotx eq 1 and iploty eq 1) then begin
  plot,x,y,xstyle=1,ystyle=1,$ 
  title=tname,xtitle='x-axis',ytitle='y-axis' 
endif else begin
  plot,x,y,/noerase,/nodata,$
  xstyle=1,ystyle=1,title=tname,xtitle='x-axis',ytitle='y-axis'
endelse

if(kmin eq kmax) then begin
 xmin = min(kkin)
 xmax = max(kkin)
endif else begin
  xmin = kmin
  xmax = kmax
endelse

mys = 255.0*(kkin - xmin)/(xmax-xmin)

if(ips eq 0) then begin		; Windows Output
  npx = long(xwid)
  npy = long(ywid)

  px0 = long((x0+lm)*npx)
  py0 = long((y0+bm)*npy)
  px1 = long((x1-rm)*npx)
  py1 = long((y1-tm)*npy)

  nx1 = long((1.0-0.4)*dx*npx)
  ny1 = long((1.0-0.5)*dy*npy)

  xtmp = size(kkin)
  nx = long(xtmp(1))
  ny = long(xtmp(2))

;  if(nx ne nx1) then begin
;   print,'Error in plotbars.pro. nx <> nx1'
;   print,nx,nx1
;   stop
;  endif

;   print,'Error in plotbars.pro. ny <> ny1'
;   print,ny,ny1
;   stop
;  endif

  if((ny ne ny1) or (nx ne nx1)) then begin
   tmp = congrid(mys,nx1,ny1)
   mys = tmp
  endif

endif else begin

  npx = 640
  npy = 512
  nx1 = long((1.0-0.4)*dx*npx)
  ny1 = long((1.0-0.5)*dy*npy)

  px0 = (x0+lm)*xwid
  py0 = (y0+bm)*ywid
  px1 = (x1-rm)*xwid
  py1 = (y1-tm)*ywid
  
  xw1 = px1-px0
  yw1 = py1-py0

endelse
 


loadct,ict

mytv,mys,px0,py0,xsize=xw1,ysize=yw1,/device

loadct,0

!p.position=[x0+lm,y0+bm,x1-rm,y1-tm]
plot,x,y,/noerase,xticklen=1.0,$
yticklen=1.0,xgridstyle=1,ygridstyle=1,xstyle=1,ystyle=1,/nodata

ndx1 = long(0.05*dx*npx)

if(iplotx eq 1) then begin
  px0 = long((x0 + 0.1*dx)*npx) 
  xa = x0+0.1*dx
  xb = x0+0.15*dx
endif else begin
  px0 = long((x1-rm+0.1*dx)*npx)
  xa = x1-rm+0.1*dx
  xb = x1-rm+0.15*dx
endelse

bar = fltarr(ndx1,ny1)
for i=0,ny1-1 do bar(*,i)=255.0*float(i)/float((ny1-1))

if(ips eq 1) then begin		; Postscript Output
  px0 = xa*xwid
  py0 = (y0+bm)*ywid
  px1 = xb*xwid
  py1 = (y1-tm)*ywid
  
  xw1 = px1-px0
  yw1 = py1-py0
endif


loadct,ict
mytv,bar,px0,py0,xsize=xw1,ysize=yw1,/device
loadct,0

y = [xmin,xmax]
x = [0,1]
!p.position=[xa,y0+bm,xb,y1-tm]
plot,x,y,/noerase,/nodata,xticks=1,$
xticklen=1e-3,xstyle=4,ystyle=1
plots,0.0,xmin
plots,1.0,xmin,/continue
plots,0.0,xmax
plots,1.0,xmax,/continue

!p.position = [x0+lm,y0+bm,x1-rm,y1-tm]

return

end


