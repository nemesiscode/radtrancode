      subroutine datarchive(aname)
C***********************************************************************
C_TITL:	DATARCHIVE
C
C_DESC:	Subroutine to insert the location of the Radtran data archive 
C	into a filename.
C
C_ARGS:	Input variables:
C	aname		CHARACTER*128	Filename of files within the
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
      INTEGER i,j,L
      CHARACTER*128 tname,arcfile
      CHARACTER aname*(*)
      LOGICAL skip,fexist
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      skip = .TRUE.

      L = LEN(ANAME)
      
      arcfile = 'datapath.arc'
      
1     format(a)
c  ** if arcfile exists then get tname from that **
      inquire(file=arcfile,exist=fexist)
      if ( fexist ) then
         open(67,file=arcfile,status='old')
          read(67,1)tname
         close(67)
            
         call remsp(tname)
         if(idiag.gt.0)then
          print*,'datarchive.f: override raddata directory to:'
          write(6,1)tname
         endif
      else

c      tname = '/Users/patirwin/radtrancode/trunk/raddata/'
<<<<<<< HEAD
c      tname ='/data/nemesis/nemesis_git/radtrancode/raddata/'
=======
      tname ='/scratch/nemesis/lnf2/radtrancode/raddata/'
>>>>>>> f9f904a4cb8c67d201fe025e2c99d910d1571631
c      tname ='/Users/glnat/CIRS/Radtran/radtrancode/raddata/'
c       tname ='/Users/irwin/gitradtran/radtrancode/raddata/'

       tname ='/network/aopp/oxpln98/plan/irwin/gitradtran/
     &radtrancode/raddata/'

      endif
      
      j = 1
      DO 10 i=1,L
        IF(skip)THEN
          IF(tname(i:i+1).eq.'/ ') THEN
            skip = .FALSE.
          ENDIF
       ELSE
         tname(i:i) = aname(j:j)
         j = j + 1
       ENDIF
         
10    CONTINUE

      aname(1:L)=tname(1:L)

      END
