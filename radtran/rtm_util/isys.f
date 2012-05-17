      INTEGER FUNCTION ISYS()
C     $Id: isys_linux.f,v 1.4 2011-01-20 14:40:06 teanby Exp $
C     ********************************************************************
C     Function to return the size of a unit of system storage for use
C     in unformatted direct access open statements etc
C     Modify to return 1, if using with OSF Compiler
C		       4, if using Fujitsu Linux compiler
C		       4, if using SunOS f90 compiler,
C		       4, if using Cray cf90 compiler (then use
C			"assign -f 90 -s unblocked -N ieee_dp" to convert
C			to 64-bit precision),
C		       4, if using Intel FORTRAN compiler v7.* for Linux,
C		       1, if using Intel FORTRAN compiler v8.0 or greater
C			for Linux. This is because Intel adopted the
C			Compaq interpretation of "numeric units".
C			Alternatively, specify "-assume byterecl" when
C			compiling.
C                      4, for GFORTRAN on Mac OSX
C		       4, if using Intel FORTRAN compile on 32-bit machines
C
C	1/1/90	SBC	Original Version.
C	3/10/94	PGJI	Updated Header.
C	28/2/00	PGJI	This version.
C	17/5/12 PGJI	Updated headers
C     ********************************************************************

      ISYS=1

      RETURN
      END
