pro readsvp,ngas,idsvp,svp

openr,2,'$HOME/radtrancode/raddata/SVP.dat'
comment = ''
loop1 : readf,2,comment
        if (strmid(comment,0,1) eq '#') then goto,loop1 
ngas = 1
reads,comment,ngas
data = fltarr(5,ngas)
readf,2,data
close,2

idsvp = lindgen(ngas)
idsvp(*) = long(data(0,*))
svp = fltarr(4,ngas)
for i=0,3 do svp(i,*)=data(i+1,*)

return

end
