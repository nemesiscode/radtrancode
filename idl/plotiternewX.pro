; ********************************************************************
; IDL program to read in and plot a Nemesis .itr file
;
; ********************************************************************

print,'Enter Nemesis run root name'
filename=' '
read,filename

spename = strcompress(filename + '.spx',/REMOVE_ALL)
apname = strcompress(filename + '.apr',/REMOVE_ALL)
prfname = strcompress(filename + '.prf',/REMOVE_ALL)
itname = strcompress(filename + '.itr',/REMOVE_ALL)

; ******* Read in measured spectrum ********************
readcirsspavX,spename,fwhm,xlat,xlon,ngeom,nav,nconv,wave,angles,spec,error

; ******* Read in .prf file to get pressure grid *******
readprfhead,prfname,npro,press

; ******* Read in .apr file to get a priori measurement vector *******
readapriori,apname,npro,nvar,varident,varparam,nx,xa,erra


phi=1.0
chisq=1.0

; ****  Open the .itr file and read in *****
openr,1,itname
 mx = 400
 nx=1
 ny=1
 niter=1
 readf,1,nx,ny,niter
 xn=fltarr(nx)
 xa=xn
 y=fltarr(ny)
 se = y
 yn=y
 yn1=y
 kk = fltarr(ny,nx) 
 tmp = fltarr(ny)
 for iter = 0,niter-1 do begin
  readf,1,chisq,phi
  readf,1,xn
  readf,1,xa
  readf,1,y
  readf,1,se
  readf,1,yn1
  readf,1,yn
  print,strcompress(string('Iteration ',iter,'  xa     xn'))
  for i=0,nx-1 do begin
   if(i le npro-1) then begin
    print,i,press(i),xa(i),exp(xa(i)),xn(i),exp(xn(i))
   endif else print,i,press(0),xa(i),exp(xa(i)),xn(i),exp(xn(i))
   
  endfor
  for i=0,nx-1 do begin
   readf,1,tmp
   kk(*,i)=tmp(*)
  endfor

  !p.multi=[0,1,2]
  plot,y(*),title='Measured and fitted y-vector',ytitle='Radiance'
  oplot,y(*)+sqrt(se(*)),linestyle=1
  oplot,y(*)-sqrt(se(*)),linestyle=1
  oplot,yn(*),linestyle=2
  oplot,yn(*),psym=1 
  oplot,yn1(*),linestyle=3
  oplot,yn1(*),psym=2
 
  plot,y(*)-yn1(*),title = 'Measured - Calculated',ytitle='Radiance'
  oplot,y(*)-yn(*),linestyle=1

  print,'Chisq, phi = ',chisq,phi
  ans=' '
  read,ans

   iprof=0
   for ivar=0,nvar-1 do begin
    itype = varident(ivar,2)
    if(itype eq 0) then iprof=iprof+1
    for i=0,nx-1 do begin
;    non-T elements are held as logs - need to convert
;     if(varident(ivar,0) ne 0)then kk(*,i)=kk(*,i)*exp(-xn(i))
    endfor
   endfor

   !p.multi=[0,iprof,1]

   istart = 0 
   for ivar=0,nvar-1 do begin
    itype = varident(ivar,2)
    print,varident(ivar,*)
    print,itype
    case itype of
     0: np=npro
     1: np=2
     2: np=1
     3: np=1
     4: np=3
     6: np=2
     8: np=3
     9: np=3
     10: np=4
     11: np=2
     888:np=long(varparam(ivar,0))
     999:np=1
     777:np=1
    endcase

    if(itype eq 0) then begin
     plot_io,xn(istart:(istart+np-1)),press,yrange=[10,0.01],$
      ytitle='Pressure (atm)',$
      title=string(varident(ivar,0),varident(ivar,1),varident(ivar,2))
     oplot,xa(istart:(istart+np-1)),press,linestyle=1
    endif

    istart=istart+np

   endfor   

   print,'plot kk matrix?'
   ans=' '
   read,ans 
   if(ans eq 'y' or ans eq 'Y')then begin 

    !p.multi=[0,1,4]
    for i=0,nx-1 do begin
     print,i,' of ',nx
;    ikeep = where(wave ge xr(0) and wave le xr(1))
;    a = max([y(ikeep),yn(ikeep)])

     plot,y,ytitle='Radiance',title='Measured and fitted radiances'
     oplot,yn,linestyle=1

     plot,yn-y,ytitle='Radiance',title='Fitted - Measured radiances'

     plot,kk(*,i),title='kk'


     plot,yn(*),$
      ytitle='Radiance',title='Radiance +/- 0.1,0.3,0.5*kk'
     oplot,yn(*)+0.1*kk(*,i),linestyle=1
     oplot,yn(*)+0.3*kk(*,i),linestyle=2
     oplot,yn(*)+0.5*kk(*,i),linestyle=3
     oplot,yn(*)-0.1*kk(*,i),linestyle=1
     oplot,yn(*)-0.3*kk(*,i),linestyle=2
     oplot,yn(*)-0.5*kk(*,i),linestyle=3
  
     read,ans

    endfor
   endif
 endfor

close,1

end
