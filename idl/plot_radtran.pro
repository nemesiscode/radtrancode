print,'Enter filename'
filename=' '
read,filename
openr,1,filename
readf,1,nplot
head=' '
readf,1,head
readf,1,xmin,xmax,ymin,ymax
readf,1,head
xname=' '
yname=' '
readf,1,xname
readf,1,yname
xname=strcompress(xname) 
yname=strcompress(yname) 
data=fltarr(2,nplot)
readf,1,data
close,1
print,'New plot (0) or overplot (1) ?'
read,ichoice
if (ichoice eq 0) then begin

	print,'xlimits and ylimits are :',xmin,xmax,ymin,ymax
        print,'Enter new values'
        read,xmin,xmax,ymin,ymax
        print,'Enter title'
        name=' '
        read,name
        plot,data(0,*),data(1,*),xrange=[xmin,xmax],yrange=[ymin,ymax],$
	xtitle=xname,ytitle=yname,title=name

endif else begin
	print,'Enter linestyle'
        read,iline
 	oplot,data(0,*),data(1,*),linestyle=iline
endelse

end
