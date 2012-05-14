      SUBROUTINE ADJUSTGRAD(dtautmpdq,amount,ngas,nlayer)
C     ********************************************************
C     Subroutine to slightly correct gradients for the case when the
C     sum of vmrs is forced to add up to 1.
C
C     Input variables
C	dtautmpdq(maxlay,maxgas+2+maxcon) Real	Rate of change of optical
C						depth with constituents at
C						each layer
C	amount(maxlay,maxgas)	Real		Layer amounts
C	ngas			Integer		Number of gases
C	nlayer			Integer		Number of layers
C
C     Output variable
C	dtautmpdq(maxlay,maxgas+2+maxcon) Real	Corrected gradients
C
C     Pat Irwin	Original	12/5/12
C
C     ********************************************************
      implicit none
      include '../includes/arrdef.f'
      integer ngas,nlayer,i,j,k
      real dtautmpdq(maxlay,maxgas+2+maxcon),tmp(maxgas)
      real amount(maxlay,maxgas)
      real totamex(maxgas),tmpam(maxgas)


      do 10 i=i,nlayer
   
C      For each gas find the sum of amounts of all other gases.
       do j=1,ngas
        totamex(j)=0.
        tmpam(j)=amount(i,j)
        do k=1,ngas
         if(k.ne.j)totamex(j)=totamex(j)+amount(i,k)
        enddo
       enddo

C      Now compute correction to gradient for each gas
       do j=1,ngas
        tmp(j)=dtautmpdq(i,j)
        do k=1,ngas
         if(k.ne.j)then
          tmp(j)=tmp(j)-tmpam(k)*dtautmpdq(i,k)/totamex(j)
         endif
        enddo
       enddo

C      Write corrected gradients back
       do j=1,ngas
        dtautmpdq(i,j)=tmp(j)
       enddo

10    continue

      return

      end
