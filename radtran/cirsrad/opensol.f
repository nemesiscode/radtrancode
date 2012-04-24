      subroutine opensol(solfile,solname)
C     ****************************************************************
C     Routine to read in solar or stellar flux reference file name
C     and fill in path to raddata directory
C  
C     Input variable
C	solfile	character*100	Name of file to read
C
C     Output variable
C 	solname	character*100	Solar/Stellar reference flux name
C
C     Pat Irwin	1/3/12	Original
C 
C     ****************************************************************
      character*100 solfile,solname
      integer ispace,ispace1

1     format(a)

      open(36,file=solfile,status='old')
        read(36,1)solname
      close(36)

C     Fill in path to the raddata directory
      call datarchive(solname)

      return

      end     
