      subroutine kk_new_sub(npoints, vi, k, vm, nm, n)
C     ***********************************************************
C     Subroutine to do Kramers-Kronig analysis and construct real part
C     of complex refractive index spectrum, give then impaginary part of the
C     refractive index spectrum and the real value at a specified wavenumber.
C 
C     Converted from IDL code written by Carly Howett
C
C     Input variables
C  	npoints		integer	Number of wavenumbers
C	vi(npoints)	real	Wavenumbers
C	k(npoints)	real	Imaginary part of RI spectrum
C	vm		real	Reference wavenumber
C	nm		real	Real part of RI at reference wavenumber
C
C     Output variable
C	n(npoints)	real	Real part of RI spectrum
C
C     Pat Irwin		24/2/14	Original
C			13/7/21 Revised and debugged
C
C     ***********************************************************
      implicit none     
      integer mpoint,i,j,l,npoints
      real pi,sum,dv
      parameter (mpoint=1000,pi=3.1415927)
      real vi(mpoint),k(mpoint),n(mpoint),vm,nm
      real va(mpoint),na(mpoint),ka(mpoint),d1,d2
      real km,alpha,beta,v,y(mpoint),delta
      logical irev   


      if(vi(1).gt.vi(npoints))then
C      Wavenumbers and spectra in reverse order
       do i=1,npoints
        va(i)=vi(npoints+1-i)
        ka(i)=k(npoints+1-i)
       enddo
       irev=.true.
      else
       do i=1,npoints
        va(i)=vi(i)
        ka(i)=k(i)
       enddo
       irev=.false.
      endif

C      do i=1,npoints
C       print*,i,va(i),ka(i)
C      enddo 

C     find km at vm
      call verint(va,ka,npoints,km,vm)

C      print*,'vm,nm,km',vm,nm,km

C     ###########################################
C     Integration
C     ###########################################     
    
                    
      do 100 i=1, npoints
    
C       sets the value of v to increment wavenumbers and presets Sum to 0
        v=va(i)
        do l=1,npoints
          y(l)=0.0
        enddo
            
        do 200 j=1, npoints
                
C          calculates the value of the 
C          potential singularity
           alpha=((va(j))**2)-(v**2)
           beta=(va(j)**2)-(vm**2)
C           print*,'alpha,beta',alpha,beta
C          if alpha isn't 0 (i.e va ne v) 
C          then K-K is evaluated for each wavenumber

           if (alpha.ne.0.and.beta.ne.0.) then
            d1 = ka(j)*va(j)-ka(i)*va(i)
            d2 = ka(j)*va(j)-km*vm
            y(j)=d1/alpha-d2/beta                        
C            print*,j,i,ka(j),ka(i),va(j),va(i),d1,d2,y(j)
           endif
                              
200     continue
            
            
C       adds up the Sum of the y values from y(1) to y(N)
C       (where N eq total no. of points)
        Sum=0.

        do l=1, npoints-1
         dv=va(l+1)-va(l)
         Sum=Sum+0.5*(y(l)+y(l+1))*dv
        enddo
                        
        na(i)=nm+(2./pi)*Sum
C        print*,i,na(i)    
100   continue      

      if(irev)then
C      reverse back output spectrum
       do i=1,npoints
        n(i)=na(npoints+1-i)  
       enddo
      else
       do i=1,npoints
        n(i)=na(i)
       enddo
      endif

      end        

        
        
    
    
    
    
