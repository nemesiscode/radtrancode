      integer function rdamform(ipfile)
C     ******************************************************************
C     Simple routine to read the header of a .prf file to see what the format
C     flag, AMFORM, is and return it.
C
C     Input variable
C	ipfile	character*100	Input filename
C
C     Output variable
C	rdamform	integer	Returned AMFORM. 0 for standard profile, 1 for
C				profiles where at each level the sum of vmrs
C				adds up to 1.
C
C     Pat Irwin	Original	12/5/12
C
C     ******************************************************************
      character*100 buffer,ipfile
      integer amform

1     format(a)
      call file(ipfile,ipfile,'prf')
      OPEN(UNIT=12,FILE=ipfile,STATUS='OLD')
C First skip the header (if any)
54    READ(12,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)AMFORM

      rdamform=AMFORM

      return

      end
