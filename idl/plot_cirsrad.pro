print,'Enter filename'
filename=' '
head=' '
xtitle=' '
ytitle=' '
xlim1=fltarr(4)
xlim2=fltarr(4)
read,filename
openr,1,filename
readf,1,n
print,'There are ',n,'spectra'
print,'Enter which one'
read,ispec
for i=1,ispec do begin
 readf,1,n
 readf,1,head
 readf,1,xlim1
 readf,1,head
 readf,1,xtitle
 readf,1,ytitle
 data1=fltarr(2,n)
 readf,1,data1
 readf,1,n
 readf,1,head
 readf,1,xlim2
 readf,1,head
 readf,1,xtitle
 readf,1,ytitle
 data2=fltarr(2,n)
 readf,1,data2
endfor
close,1

!p.multi=[0,1,2]

plot,data1(0,*),data1(1,*),xrange=xlim1(0:1),yrange=[0,xlim1(3)],$
xtitle='Wavenumber (cm-1)',title='Calculation spectrum'

plot,data2(0,*),data2(1,*),xrange=xlim2(0:1),yrange=[0,xlim2(3)],$
xtitle='Wavenumber (cm-1)',title='Convolved spectrum'


end
