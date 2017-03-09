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
readprfhead4,prfname,npro,press

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
  print,'i,xa(i),exp(xa(i)),xn(i),exp(xn(i)'
  for i=0,nx-1 do begin
    print,i,xa(i),exp(xa(i)),xn(i),exp(xn(i))
;   if(i le npro-1) then begin
;    print,i,press(i),xa(i),exp(xa(i)),xn(i),exp(xn(i))
;   endif else print,i,press(0),xa(i),exp(xa(i)),xn(i),exp(xn(i))
   
  endfor
  for i=0,nx-1 do begin
   readf,1,tmp
   kk(*,i)=tmp(*)
  endfor

  window,0,xsize=1800,ysize=800
  !p.multi=[0,1,3]
  a = min([y(*)-sqrt(se(*)),yn(*)])
  b = max([y(*)+sqrt(se(*)),yn(*)])

;  plot,wave,y(*),title='Measured and fitted y-vector',ytitle='Radiance',$
;    yrange=[a,b]
;  oplot,wave,y(*)+sqrt(se(*)),linestyle=1
;  oplot,wave,y(*)-sqrt(se(*)),linestyle=1
;  oplot,wave,yn(*),linestyle=2
;  oplot,wave,yn(*),psym=1 
;  oplot,wave,yn1(*),linestyle=3
;  oplot,wave,yn1(*),psym=2

  plot,y(*),title='Measured and fitted y-vector',ytitle='Radiance',$
    yrange=[a,b]
  oplot,y(*)+sqrt(se(*)),linestyle=1
  oplot,y(*)-sqrt(se(*)),linestyle=1
  oplot,yn(*),linestyle=2
;  oplot,yn(*),psym=1 
  oplot,yn1(*),linestyle=3
;  oplot,yn1(*),psym=2


  a = min(yn(*))
  b = max(yn(*))

;  plot_io,wave,y(*),title='Measured and fitted y-vector',ytitle='Radiance',$
;    yrange=[a,b]
;  oplot,wave,y(*)+sqrt(se(*)),linestyle=1
;  oplot,wave,y(*)-sqrt(se(*)),linestyle=1
;  oplot,wave,yn(*),linestyle=2
;  oplot,wave,yn(*),psym=1 
;  oplot,wave,yn1(*),linestyle=3
;  oplot,wave,yn1(*),psym=2

  plot_io,y(*),title='Measured and fitted y-vector',ytitle='Radiance',$
    yrange=[a,b]
  oplot,y(*)+sqrt(se(*)),linestyle=1
  oplot,y(*)-sqrt(se(*)),linestyle=1
  oplot,yn(*),linestyle=2
;  oplot,yn(*),psym=1 
  oplot,yn1(*),linestyle=3
;  oplot,yn1(*),psym=2

 
;  plot,wave,y(*)-yn1(*),title = 'Measured - Calculated',ytitle='Radiance'
;  oplot,wave,y(*)-yn(*),linestyle=1

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
     12: np=3
     13: np=3
     444:np=2+long(varparam(ivar,0))
     555:np=1
     666:np=1
     888:np=long(varparam(ivar,0))
     887:np=long(varparam(ivar,0))
     889:np=1
     999:np=1
     222:np=8
     223:np=9
     224:np=9
     225:np=11
     777:np=1
    endcase

    if(itype eq 0) then begin
     if(varident(ivar,0) ne 0) then begin
      plot_oo,exp(xn(istart:(istart+np-1))),press,yrange=[10,0.01],$
       ytitle='Pressure (atm)',$
       title=string(varident(ivar,0),varident(ivar,1),varident(ivar,2))
       oplot,exp(xa(istart:(istart+np-1))),press,linestyle=1
     endif else begin
      plot_io,xn(istart:(istart+np-1)),press,yrange=[10,0.01],$
       ytitle='Pressure (atm)',$
       title=string(varident(ivar,0),varident(ivar,1),varident(ivar,2))
       oplot,xa(istart:(istart+np-1)),press,linestyle=1
     endelse
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

     plot,y,ytitle='Radiance',title='Measured and fitted radiances',yrange=[a,b]
     oplot,yn,linestyle=1

     plot,yn-y,ytitle='Radiance',title='Fitted - Measured radiances'

     plot,kk(*,i),title='kk'


     plot,yn(*),$
      ytitle='Radiance',title='Radiance +/- 0.1,1,10,100*kk',yrange=[a,b]
     oplot,yn(*)+0.1*kk(*,i),linestyle=1
     oplot,yn(*)+kk(*,i),linestyle=2
     oplot,yn(*)+10*kk(*,i),linestyle=3
     oplot,yn(*)+100*kk(*,i),linestyle=4
     oplot,yn(*)-0.1*kk(*,i),linestyle=1
     oplot,yn(*)-kk(*,i),linestyle=2
     oplot,yn(*)-10*kk(*,i),linestyle=3
     oplot,yn(*)-100*kk(*,i),linestyle=4
  
     read,ans

    endfor
   endif
 endfor

close,1

end
