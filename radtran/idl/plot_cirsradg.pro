print,'Enter filename'
filename=' '
gname=strarr(30)
tname=' '
npath=1
nconv=1
nv=1
read,filename
openr,1,filename
readf,1,npath,nconv,nv
for iv=0,nv-1 do begin
 readf,1,tname
 gname(iv)=tname
endfor
print,'There are ',npath,' calculations'
print,'Enter which one'
read,ipath
for i=1,ipath do begin
 n=1
 readf,1,n
 spec=fltarr(2,nconv)
 readf,1,spec
 grads = fltarr(nconv,nv)
 readf,1,grads
endfor
close,1

!p.multi=[0,1,2]

for iv=0,nv-1 do begin
 plot,spec(0,*),spec(1,*),xtitle='Wavenumber (cm-1)',$
 title='Calculated spectrum'

 plot,spec(0,*),grads(*,iv),xtitle='Wavenumber (cm-1)',$
  title='Gradient : ' + gname(iv)

 ans=' '
 read,ans

endfor

end
