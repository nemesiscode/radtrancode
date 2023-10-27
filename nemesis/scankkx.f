      subroutine scankkx(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  varparam,kkx,stx,nxx)
C     *****************************************************************
C     Subroutine to strip elements for current retrieval 
C     parameters from Jacobian calculated for previously retrieved 
C     parameters. Part of LIN=3 option where previously
C     retrieved parameters are used to set apriori for same parameters in
C     current retrieval and used to fix profiles in .prf files for other  
C     parameters not being retrieved this time round.
C
C     Input variables
C      nvarx integer Number of parameters previously retrieved
C                    (and listed in .pre file)
C      varidentx(mvar,3) integer parameter specifiers
C      varparamx(mvar,mparam) real parameter details 
C      nprox integer Number of vertical levels in previous .prf files
C      nvar integer Number of parameters currently being retrieved
C      varident(mvar,3) integer parameter specifiers
C
C     Output variables
C      kkx(my,mx) real Jacobian calculated for previously retrieved vector
C                      where elements currently being retrieved have been
C                      stripped out.
C      stx(mx,mx) real State covariance matrix where elements currently 
C                      being retrieved have been stripped out.

C      nxx integer Number of elements in stripped retrieval vector
C
C     Pat Irwin   7/6/06
C
C     *****************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer nprox,nvarx,varidentx(mvar,3),nxx,ivar,ivarx,npvar
      real varparamx(mvar,mparam)
      integer npro,nvar,varident(mvar,3),ioff,ikeep,np,icopy,i,j
      real varparam(mvar,mparam)
      real kkx(my,mx),kk1(my,mx),stx(mx,mx),st1(mx,mx)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C     **************************** CODE ********************************

      ioff = 0
      ikeep = 0

      do i=1,mx
       do j=1,my
        kk1(j,i)=0.0
       enddo
       do j=1,mx
        st1(j,i)=0.0
       enddo
      enddo

      do ivarx = 1,nvarx
        np=-1
        if(varidentx(ivarx,1).le.100)then
          np = npvar(varidentx(ivarx,3),nprox,varparamx(ivarx,1))
        else
         if(varidentx(ivarx,1).eq.102)np = 1
         if(varidentx(ivarx,1).eq.555)np = 1
         if(varidentx(ivarx,1).eq.333)np = 1
         if(varidentx(ivarx,1).eq.222)np = 8
         if(varidentx(ivarx,1).eq.223)np = 9
         if(varidentx(ivarx,1).eq.224)np = 9
         if(varidentx(ivarx,1).eq.225)np = 11
         if(varidentx(ivarx,1).eq.226)np = 8
         if(varidentx(ivarx,1).eq.227)np = 7
         if(varidentx(ivarx,1).eq.999)np = 1
         if(varidentx(ivarx,1).eq.888)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,1).eq.887)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,1).eq.444)then
          if(varparamx(ivarx,2).gt.0.0)then
           np = 2+int(varparamx(ivarx,1))
          else
           np = 3
          endif
         endif
         if(varidentx(ivarx,1).eq.446)then
          if(varparamx(ivarx,2).gt.0.0)then
           np = 3+2*int(varparamx(ivarx,1))
          else
           np = 5
          endif
         endif
         if(varidentx(ivarx,1).eq.445)then
          np = 3+(2*int(varparamx(ivarx,1)))
         endif
         if(varidentx(ivarx,1).eq.777)np = 1
         if(varidentx(ivarx,1).eq.666)np = 1
        endif
        icopy = 1
        do ivar = 1,nvar

          if(varidentx(ivarx,1).eq.varident(ivar,1))then
           if(varidentx(ivarx,2).eq.varident(ivar,2))then            
            icopy=0
             if(varidentx(ivarx,3).eq.28)then
              if(varparamx(ivarx,1).ne.varparam(ivar,1))then
               icopy=1
              endif
             endif
           endif
          endif
        enddo


        if(icopy.gt.0)then
         do i=1,np
          do j=1,my
           kk1(j,ikeep+i)=kkx(j,ioff+i)
          enddo
          do j=1,mx
           st1(j,ikeep+i)=stx(j,ioff+i)
          enddo
         enddo
         ikeep=ikeep+np
        endif

        ioff = ioff+np

      enddo

      nxx = ikeep

      do i=1,mx
       do j=1,my
        kkx(j,i)=kk1(j,i)
       enddo
       do j=1,mx
        stx(j,i)=st1(j,i)
       enddo
      enddo

      return

      end
