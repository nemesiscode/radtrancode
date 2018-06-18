      SUBROUTINE ADD_LINE_BUFF(BUFFER,NLIN)
      INTEGER NLIN
      CHARACTER*256 BUFFER
      INCLUDE '../includes/lincomseq.f'

      IF(NLIN.GT.MAXLINSEQ)THEN
       PRINT*,'Error in add_line_buff.f NLIN > MAXLINSEQ'
       PRINT*,NLIN,MAXLINSEQ
       STOP
      ENDIF

      SEQBUF(NLIN)=BUFFER

      RETURN

      END
