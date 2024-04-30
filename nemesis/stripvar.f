      subroutine stripvar(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  varparam,nxx,xnx)
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
      real varparamx(mvar,mparam),xnx(mx),varparam(mvar,mparam)
      integer nvar,varident(mvar,3),ikeep,np,icopy,i,j,ioff,nlong
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C     **************************** CODE ********************************

      
      ikeep = 0
      ioff = 0
      ivar=0
      do ivarx = 1,nvarx
        np=1
        if(varidentx(ivarx,1).le.100)then
         if(varidentx(ivarx,3).eq.0)np = nprox
         if(varidentx(ivarx,3).eq.1)np = 2
         if(varidentx(ivarx,3).eq.2)np = 1
         if(varidentx(ivarx,3).eq.3)np = 1
         if(varidentx(ivarx,3).eq.4)np = 3
         if(varidentx(ivarx,3).eq.7)np = 2
         if(varidentx(ivarx,3).eq.8)np = 3
         if(varidentx(ivarx,3).eq.9)np = 3
         if(varidentx(ivarx,3).eq.7)np = 2
         if(varidentx(ivarx,3).eq.10)np = 4
         if(varidentx(ivarx,3).eq.11)np = 2
         if(varidentx(ivarx,3).eq.12)np = 3
         if(varidentx(ivarx,3).eq.13)np = 3
         if(varidentx(ivarx,3).eq.14)np = 3
         if(varidentx(ivarx,3).eq.15)np = 3
         if(varidentx(ivarx,3).eq.16)np = 4
         if(varidentx(ivarx,3).eq.17)np = 2
         if(varidentx(ivarx,3).eq.18)np = 2
         if(varidentx(ivarx,3).eq.19)np = 4
         if(varidentx(ivarx,3).eq.20)np = 2
         if(varidentx(ivarx,3).eq.21)np = 2
         if(varidentx(ivarx,3).eq.22)np = 5
         if(varidentx(ivarx,3).eq.23)np = 4
         if(varidentx(ivarx,3).eq.24)np = 3
         if(varidentx(ivarx,3).eq.25)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,3).eq.26)np = 4
         if(varidentx(ivarx,3).eq.27)np = 3
         if(varidentx(ivarx,3).eq.28)np = 1
         if(varidentx(ivarx,3).eq.29)np = int(varparamx(ivarx,1))*nprox
         if(varidentx(ivarx,3).eq.30)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,3).eq.31)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,3).eq.32)np = 3
         if(varidentx(ivarx,3).eq.33)then
           NLONG = INT(VARPARAMX(IVARX,1)/VARPARAMX(IVARX,2)+0.1)
           np = int(varparamx(ivarx,1)) + NLONG
         endif
         if(varidentx(ivarx,3).eq.34)np = 2
         if(varidentx(ivarx,3).eq.35)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,3).eq.36)np = int(varparamx(ivarx,1))
         if(varidentx(ivarx,3).eq.37)np = 1
         if(varidentx(ivarx,3).eq.38)np = 1
         if(varidentx(ivarx,3).eq.39)np = 1
         if(varidentx(ivarx,3).eq.40)np = 2
         if(varidentx(ivarx,3).eq.41)np = 5
         if(varidentx(ivarx,3).eq.42)np = 3
         if(varidentx(ivarx,3).eq.43)np = 5
         if(varidentx(ivarx,3).eq.44)np = 5*int(varparamx(ivarx,1))
         if(varidentx(ivarx,3).eq.45)np = 3
         if(varidentx(ivarx,3).eq.46)np = 6
         if(varidentx(ivarx,3).eq.47)np = 3
         if(varidentx(ivarx,3).eq.48)np = 4
         if(varidentx(ivarx,3).eq.49)np = 1
         if(varidentx(ivarx,3).eq.50)np = 3
         if(varidentx(ivarx,3).eq.51)np = 4
         if(varidentx(ivarx,3).eq.52)np = 4
        else
         if(varidentx(ivarx,1).eq.102)np = 1
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
             if(varidentx(ivarx,3).eq.28)then
              if(varparamx(ivarx,1).ne.varparam(i,1))then
               icopy=1
              endif
             endif
           endif
          endif
        enddo

        if(idiag.gt.0)print*,'ivarx',ivarx
        if(idiag.gt.0)print*,'icopy',icopy

        if(icopy.gt.0)then
          ivar = ivar+1
         
          do j=1,3
           varidentx(ivar,j)=varidentx(ivarx,j)
          enddo
          do j=1,mparam
           varparamx(ivar,j)=varparamx(ivarx,j)
          enddo
          do j=1,np
           xnx(ikeep+j)=xnx(ioff+j)
          enddo

          ikeep=ikeep+np

        endif

        ioff=ioff+np

      enddo

      nxx = ikeep
      nvarx = ivar

      return

      end
