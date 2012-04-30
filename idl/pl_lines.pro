; *******************************************************************
; IDL routine to plot the ASCII output of the new RADTRAN program
; Pl_alines
;
; Pat Irwin	13/1/00
; *******************************************************************

filename = ' '
print,'Enter filename'
read,filename
openr,1,filename
xname=' '
gname=' '
dname=' '
head=' '
readf,1,xname
readf,1,gname
readf,1,dname
readf,1,head
nline=1
readf,1,nline

lines=fltarr(2,nline)
readf,1,lines

close,1
tname = strcompress(gname+dname)

iplot=0
print,'New plot(0) or overplot (1)?'
read,iplot

if(iplot eq 0) then begin
 a = min(lines(1,*))
 b = max(lines(1,*))
 print,a,b
 plot,lines(0,*),lines(1,*),xtitle=xname,ytitle='Line strength',$
 title=tname,/nodata,yrange=[a,b],ystyle=1
endif

print,'Enter linestyle'
read,iline

for i=0,nline-1 do begin
 plots,lines(0,i),a,linestyle=iline
 plots,lines(0,i),lines(1,i),/continue,linestyle=iline
endfor

end
