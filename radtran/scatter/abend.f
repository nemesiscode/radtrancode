	subroutine abend( msg)
C
C  print error message and exit
C
C_HIST:	???????	???	Original version.

	CHARACTER*(*) msg

	l = LEN(msg)
	IF (l.GT.0) PRINT*, msg(:l)
	PRINT*,' ABEND.f :: Stopping program.'
	STOP
cc	call exit

	END
