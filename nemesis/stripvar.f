      subroutine stripvar(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  nxx,xnx)
C     *****************************************************************
C     Subroutine to strip current retrieval parameters from previous 
C     retrieved parameters. Part of LIN=3 option where previously
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
C      nxx integer Number of elements in previous retrieval vector
C      xnx(mx) real Previously retrieved vector
C 
C     Output variables
C      nvarx integer Number of parameters previously retrieved
C                    and not beginh retrieved now.
C      varidentx(mvar,3) integer parameter specifiers
C      varparam(mvar,mparam) real parameter details 
C      nxx integer Number of elements in previous retrieval vector not
C                  beign retrieved now
C      xnx(mx) real Previously retrieved vector with currently retrieved
C                   elements stripped out.
C
C     Pat Irwin   7/6/06
C
C     *****************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer nprox,nvarx,varidentx(mvar,3),nxx,ivar,ivarx
      real varparamx(mvar,mparam),xnx(mx)
      integer nvar,varident(mvar,3),ikeep,np,icopy,i,j,ioff

C     **************************** CODE ********************************

      
      ikeep = 0
      ioff = 0
      ivar=0
      do ivarx = 1,nvarx
        np=-1
        if(varidentx(ivarx,1).le.100)then
         if(varidentx(ivarx,3).eq.0)np = nprox
         if(varidentx(ivarx,3).eq.1)np = 2
         if(varidentx(ivarx,3).eq.2)np = 1
         if(varidentx(ivarx,3).eq.3)np = 1
         if(varidentx(ivarx,3).eq.4)np = 3
         if(varidentx(ivarx,3).eq.8)np = 3
         if(varidentx(ivarx,3).eq.9)np = 3
         if(varidentx(ivarx,3).eq.7)np = 2
         if(varidentx(ivarx,3).eq.17)np = 2
         if(varidentx(ivarx,3).eq.18)np = 2
         if(varidentx(ivarx,3).eq.19)np = 4
         if(varidentx(ivarx,3).eq.20)np = 2
         if(varidentx(ivarx,3).eq.23)np = 4
         if(varidentx(ivarx,3).eq.24)np = 3
         if(varidentx(ivarx,3).eq.25)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,3).eq.26)np = 4
         if(varidentx(ivarx,3).eq.27)np = 3
        else
         if(varidentx(ivarx,1).eq.555)np = 1
         if(varidentx(ivarx,1).eq.333)np = 1
         if(varidentx(ivarx,1).eq.999)np = 1
         if(varidentx(ivarx,1).eq.888)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,1).eq.777)np = 1
         if(varidentx(ivarx,1).eq.666)np = 1
        endif
        icopy = 1
        do i = 1,nvar
          if(varidentx(ivarx,1).eq.varident(i,1))then
           if(varidentx(ivarx,2).eq.varident(i,2))then
            icopy=0
           endif
          endif
        enddo

        if(icopy.gt.0)then
          ivar = ivar+1
         
          do i=1,3
           varidentx(ivar,i)=varidentx(ivarx,i)
          enddo
          do i=1,mparam
           varparamx(ivar,i)=varparamx(ivarx,i)
          enddo
          do i=1,np
           xnx(ikeep+i)=xnx(ioff+i)
          enddo

          ikeep=ikeep+np

        endif

        ioff=ioff+np

      enddo

      nxx = ikeep
      nvarx = ivar

      return

      end
