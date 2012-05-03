      subroutine datarchive(aname)
C***********************************************************************
C_TITL:	DATARCHIVE
C
C_DESC:	Subroutine to insert the location of the Radtran data archive 
C	into a filename.
C
C_ARGS:	Input variables:
C	aname		CHARACTER*100	Filename of files within the
C					raddata directory.
C
C_FILE:	No files opened.
C
C_CALL:	No calls.
C
C_HIST:	14/8/00	PGJI
C	3aug05	NT	look for arcfile to tell progs where to look for
C				raddata archives.
C***********************************************************************

      IMPLICIT NONE

      INTEGER i,j
      CHARACTER*100 aname,tname,arcfile
      LOGICAL skip,fexist

      skip = .TRUE.
      
      arcfile = 'datapath.arc'
      
c  ** if arcfile exists then get tname from that **
      inquire(file=arcfile,exist=fexist)
      if ( fexist ) then
         open(67,file=arcfile,status='old')
         read(67,*) tname
         close(67)
      else
      tname ='/home/oxpln98/plan/irwin/radtrancode/trunks/
     &raddata/'
      endif
      
c      print*,'tname=',tname

      j = 1
      DO 10 i=1,80
        IF(skip)THEN
          IF(tname(i:i+1).eq.'/ ') THEN
            skip = .FALSE.
          ENDIF
       ELSE
         tname(i:i) = aname(j:j)
         j = j + 1
       ENDIF
         
10    CONTINUE

      aname(1:80)=tname(1:80)

      END
