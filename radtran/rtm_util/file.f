      SUBROUTINE FILE(NAME1,NAME2,EXT)
C     $Id: file.f,v 1.3 2007-07-18 11:37:41 irwin Exp $
C     ********************************************************************
C     Forces correct VMS style file extension for a filename. i.e. assumes 
C     a <4 character extension after a dot separator. 
C
C     1/1/90    Original Version:       SBC
C     3/10/94   Updated Header          PGJI
C
C     ********************************************************************
      INTEGER I,L,LE
      CHARACTER*(*) NAME1,NAME2,EXT

C     first copy the file name to the second array and remove any leading
C     spaces
      NAME2=NAME1
      CALL REMSP(NAME2)
C     now stepping back from the end of the name through the last four
C     characters
      LE=LEN(NAME2)
      IF(LE.LT.1)RETURN
      DO 40 L=LE,1,-1
      IF(NAME2(L:L).NE.' ')THEN
        LE=L
        GOTO 50
        END IF
40    CONTINUE
      RETURN
50    CONTINUE
      L=LE-3
      IF(L.LT.1)L=1
      DO 10 I=LE,L,-1
C     if a directory tree delimeter is found just add the extension to the end
C     of the name
      IF(NAME2(I:I).EQ.'/'.OR.NAME2(I:I).EQ.']')GOTO 30
C     if a dot is found overwrite any existing extension
      IF(NAME2(I:I).EQ.'.')THEN
        NAME2(I+1:I+3)=EXT
        GOTO 20
        END IF
10    CONTINUE
C     if no dot is found in last four characters add the extension to the end
C     of the filename
30    NAME2(LE+1:LE+1)='.'
      NAME2(LE+2:LE+4)=EXT

      CALL REMSP(NAME2)
     
20    RETURN
      END 
