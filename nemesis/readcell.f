      subroutine readcell(cellfile)
C     *****************************************************************
C     Subroutine to read in the properties of an SCR cell file to be 
C     added to subsequent .pat file and calculations. Data read into a 
C     common block.
C     The cell is assumed to be a Selective Chopper Radiometer (SCR) cell.
C
C     Pat Irwin		1/4/14
C
C     *****************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'

      character*100 cellfile,head
      integer cellngas,cellid(maxgas),celliso(maxgas),icread,i
      real cellength,cellpress,celltemp,cellvmr(maxgas)

      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp

      open(12,file=cellfile,status='old')
      read(12,1)head
1     format(a)
      read(12,1)head
      read(head(6:),*)cellngas
      do 10 i=1,cellngas
       read(12,*)cellid(i),celliso(i),cellvmr(i)
10    continue
      read(12,1)head
      read(head(7:),*)cellength
      read(12,1)head
      read(12,*)cellpress,celltemp

      close(12)

      return

      end
