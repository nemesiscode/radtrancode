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
      INTEGER i,j,L,evstatus,evlen
      CHARACTER*128 tname,arcfile,radrepo
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
c       If there is not `datapath.arc` fle, then assume there is an
c       environment variable called "RADREPO" that points to the
c       location of the NEMESIS repository

c       Setting `evlen=0` copies the whole environment variable
        evlen=0
        CALL GET_ENVIRONMENT_VARIABLE("RADREPO", radrepo, evlen,        &
     &evstatus, .TRUE.)

c       IF `evstatus` is zero then everything worked, otherwise
c       something went wrong
        if ( evstatus .eq. 0 ) then
          tname = TRIM(ADJUSTL(radrepo)) // '/raddata/'
        else
          if (evstatus .eq. -1 ) then
            write(*,*) 'datarchive.f:'
            write(*,*) 'Envionment variable RADREPO too long (more than &
     &128 characters).'
          endif
          if ( evstatus .eq. 1 ) then
            write(*,*) 'datarchive.f:'
            write(*,*) 'Evironment variable RADREPO does not exist.'
          endif
          if ( evstatus .eq. 2 ) then
            write(*,*) 'datarchive.f:'
            write(*,*) 'Processor does not support environment variables&
     &, cannot read RADREPO.'
          endif
          write(*,*) 'datarchive.f:'
          write(*,*) 'Could not automatically find location for ".../rad&
     &data", falling back on hard-coded value. If this fails, alter the &
     &hard coded value for "tname" in the file .../radtran/rtm_util/data&
     &rchive.f'


c         --------------------------------------------------------------
c         If all esle fails, fall back on hard-coded value. It should be
c         set to '.../raddata', where '...' is the path to the NEMESIS
c         repository on the machine NEMESIS will be running on.
c         --------------------------------------------------------------

c          tname ='/Users/patirwin/radtrancode/trunk/raddata/'
c          tname ='/data/nemesis/nemesis_git/radtrancode/raddata/'
          tname ='/Users/glnat/CIRS/Radtran/radtrancode/raddata/'
c          tname ='/Users/irwin/gitradtran/radtrancode/raddata/'
c          tname ='/network/aopp/oxpln98/plan/irwin/gitradtran/&radtranco&
c     &de/raddata/'
        endif
      endif
      
c     Concatenate 'tname' and 'aname', this should result in 
c     <path to .../raddata> and <some file in raddata we want>
c     being joined together.
      j = 1
      DO 10 i=1,L
        IF(skip)THEN
          IF(tname(i:i+1).eq.'/ ') THEN
c           This could fail if a directory name begins with a space
c           unlikely, but possible.
            skip = .FALSE.
          ENDIF
       ELSE
         tname(i:i) = aname(j:j)
         j = j + 1
       ENDIF
         
10    CONTINUE

      aname(1:L)=tname(1:L)

      END
