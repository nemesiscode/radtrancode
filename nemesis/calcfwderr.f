      subroutine calcfwderr(nx,ny,kkx,sxx,sf)
C     *************************************************************
C     Routine to calculate effective forward-modelling error of previous
C     retrievals
C
C     Pat Irwin		21/2/05
C	
C     *************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      integer nx,ny,i,j
      real kkx(my,mx),sxx(mx,mx),sf(my,my)
      double precision kk(my,my),kt(my,my),a(my,my),m(my,my)
      integer kk1,kk2,kt1,kt2,m1,m2,a1,a2

      if (nx.eq.0) then
        do i=1,ny
         do j=1,ny
          sf(i,j)=0.0
         enddo
        enddo
        return
      endif

      do i=1,ny
       do j=1,nx
        kk(i,j) = dble(kkx(i,j))
        kt(j,i) = dble(kkx(i,j))
       end do
      end do
      kk1=ny
      kk2=nx
      kt1=nx
      kt2=ny

C     Load a with sx_in

      do i=1,nx
         do j=1,nx
          a(i,j)=dble(sxx(i,j))
         end do
      end do

C     Multiply sx*kt, put answer in m
      call dmult_mat(my,a,nx,nx,kt,kt1,kt2,m,m1,m2)

C     Multiply kk*m, put answer in a
      call dmult_mat(my,kk,kk1,kk2,m,m1,m2,a,a1,a2)

      do i=1,ny
       do j=1,ny
        sf(i,j)=sngl(a(i,j))
       enddo
      enddo

      return

      end
