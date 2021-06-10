function calc_int,k
; performs numerical integration of (cos(thet))^2*k dthet, from 0 to pi/2
;
nthet=101
thet = 0.5*!pi*findgen(nthet)/float(nthet-1)

cthet = cos(thet)
cthet(nthet-1)=0.0

dthet = 0.5*!pi/float(nthet-1)

y = cthet^(2.0*k)

sum=0.0
for i=0,nthet-2 do sum = sum+0.5*(y(i)+y(i+1))*dthet

return,sum

end

