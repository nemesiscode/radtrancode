print,'Enter filename'
filename=' '
read,filename
openr,1,filename
readf,1,n
print,'There are ',n,'spectra'
print,'Enter which one'
read,ispec
for i=1,ispec do begin
 readf,1,n
 data1=fltarr(2,n)
 readf,1,data1
 readf,1,n
 data2=fltarr(2,n)
 readf,1,data2
endfor
close,1
end
