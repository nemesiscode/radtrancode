      subroutine write_covariance(runname,npro,nvar,varident,
     1 varparam,nx,ny,sa,sm,sn,st,se,aa,dd,kk)
C     **************************************************************
C     Subroutine to dump diagnostic matrices to output file for
C     inspection and assessment.
C 
C     Pat Irwin	Original	9/9/04
C
C     **************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      integer npro,nvar,varident(mvar,3),nx,ny,i,j
      real varparam(mvar,mparam),sa(mx,mx),sm(mx,mx),sn(mx,mx)
      real st(mx,mx),kk(my,mx),se(my)
      double precision aa(mx,mx),dd(mx,my)
      character*100 runname

      call file(runname,runname,'cov')

      open(33,file=runname,status='unknown')
       write(33,*)npro,nvar
       do i=1,nvar
        write(33,*)(varident(i,j),j=1,3)
        write(33,*)(varparam(i,j),j=1,5)
       enddo
       write(33,*)nx,ny

       do i=1,nx   
        write(33,*)(sa(i,j),j=1,nx)
        write(33,*)(sm(i,j),j=1,nx)
        write(33,*)(sn(i,j),j=1,nx)
        write(33,*)(st(i,j),j=1,nx)
       enddo
        
       do i=1,nx
        write(33,*)(sngl(aa(i,j)),j=1,nx)
       enddo
       
       do i=1,ny
        write(33,*)(sngl(dd(j,i)),j=1,nx)
       enddo
        
       do i=1,nx
        write(33,*)(kk(j,i),j=1,ny)
       enddo
        
       write(33,*)(se(i),i=1,ny)

      close(33)
     
      return

      end
