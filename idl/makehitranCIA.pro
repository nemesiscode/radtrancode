; 1............H2-H2 (ortho:para = 1:1 `equilibrium')
; 2............H2-He                        "
; 3............H2-H2 (ortho:para = 3:1 `normal')
; 4............H2-He                        "
; 5............H2-N2
; 6............N2-CH4
; 7............N2-N2
; 8............CH4-CH4
; 9............H2-CH4

ntempk=25
nwavek=1501
dnu=10.

vv=dnu*findgen(nwavek)
tempk=50+15*findgen(ntempk)
knew=fltarr(9,ntempk,nwavek)
knew(*)=0.

dirname='/home/valhalla/plan/sinclair/hit12/CIA/'

list1=['Main-Folder/H2-H2/H2-H2_2011.cia','Main-Folder/H2-He/H2-He_2011.cia',$
       'Main-Folder/H2-H2/H2-H2_2011.cia','Main-Folder/H2-He/H2-He_2011.cia',$
       'Main-Folder/N2-H2/N2-H2_2011.cia','Alternate-Folder/N2-CH4/N2-CH4_2011.cia',$
       'Main-Folder/N2-N2/N2-N2_2011.cia','Alternate-Folder/CH4-CH4/CH4-CH4_2011.cia',$
       'Main-Folder/H2-CH4/H2-CH4_norm_2011.cia']

ntemp1=[113,334,$
        113,334,$
        10,10,$
        10,10,$
        10]


for ipair=0,8 do begin
 filename=dirname+list1(ipair)
 ntemp=ntemp1(ipair)
 readhitranCIA,filename,ntemp,nwave,temp,wave,kkin
 print,list1(ipair)
 print,ntemp
 print,temp
 print,nwave
 print,wave
 print,'Converting...'
 convertk,ntemp,nwave,temp,wave,kkin,ntempk,nwavek,tempk,vv,kkout
 for i=0,ntempk-1 do knew(ipair,i,*)=kkout(i,*)
endfor

temps=double(tempk)
openw,1,'makehitranCIA.tab',/f77_unformatted
 writeu,1,temps
 writeu,1,knew
close,1

dnu=1.
vv=dnu*findgen(nwavek)
knew(*)=0.

for ipair=0,8 do begin
 filename=dirname+list1(ipair)
 ntemp=ntemp1(ipair)
 readhitranCIA,filename,ntemp,nwave,temp,wave,kkin
 print,list1(ipair)
 print,ntemp
 print,temp
 print,nwave
 print,wave
 convertk,ntemp,nwave,temp,wave,kkin,ntempk,nwavek,tempk,vv,kkout
 for i=0,ntempk-1 do knew(ipair,i,*)=kkout(i,*)
endfor

openw,1,'makehitranCIAiso.tab',/f77_unformatted
 writeu,1,temps
 writeu,1,knew
close,1


end

