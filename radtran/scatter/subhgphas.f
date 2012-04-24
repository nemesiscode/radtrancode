       subroutine subhgphas(nphase,theta,x,cphase,kk)
       integer max_thet,nphase
       parameter(max_thet=100)
       real theta(max_thet),x(3),cphase(max_thet),kk(max_thet,3)
       real pi,f,g1,g2,henyey,tphase(max_thet),alpha(max_thet)
       real xt(3)
       parameter (pi=3.1415927)

C      print*,x
       f = x(1)
       g1 = x(2)
       g2 = x(3)

       do 10 i=1,nphase
        alpha(i) = cos(theta(i)*pi/180.0)
        calpha = alpha(i)
        cphase(i) = henyey(calpha,f,g1,g2)
10     continue

       xt(1)=x(1)
       xt(2)=x(2)
       xt(3)=x(3)

       do 20 j=1,3
        dx = 0.01
        xt(j) = x(j)+dx
        if(j.eq.1)then
         if(xt(j).gt.0.99)xt(j)=x(j)-dx
        else if (j.eq.2)then
         if(xt(j).gt.0.98)xt(j)=x(j)-dx
        endif

        dx = xt(j)-x(j)

        f = xt(1)
        g1 = xt(2)
        g2 = xt(3)
C        print*,xt,dx
        do 15 i=1,nphase
         tphase(i) = henyey(alpha(i),f,g1,g2)
15      continue
	do 16 i=1,nphase
         kk(i,j)=(tphase(i)-cphase(i))/dx
16	continue
       
        xt(j)=x(j)
        
20     continue

       return

       end
