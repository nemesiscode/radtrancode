      subroutine forwarderr(ipfile,ngeom,nconv,vconv,woff,rerr)
C     Id:
C     *****************************************************************
C     Subroutine which returns the forward modelling error read in
C     from an external file
C
C     Input variables
C	ipfile 	character*100 forward modelling error filename
C	ngeom	integer	Number of observing geometries
C	nconv(mgeom) integer Number of convolution wavelengths
C	vconv(mgeom,mconv) real	Convolution wavelengths
C	woff	real	Any wavenumber offset applied to spectrum
C
C     Output variables
C	rerr(mgeom,mconv) real	%radiance forward-modelling errors
C 
C     Pat Irwin	4/5/01		Original
C     Pat Irwin 17/10/03	Tidied for Nemesis
C     Pat Irwin 11/8/04		Overhauled to read in errors from
C				  external file
C
C     *****************************************************************

      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer nconv(mgeom),np,i,j,ngeom,igeom
      real vconv(mgeom,mconv),rerr(mgeom,mconv),v1(mconv)
      real r1(mconv),woff,x,y
      character*100 IPFILE

      OPEN(12,FILE=IPFILE,STATUS='OLD')
      READ(12,*)NP
      DO I=1,NP
        READ(12,*)V1(I),R1(I)
        V1(I)=V1(I)+WOFF
      ENDDO
      CLOSE(12)

      DO IGEOM=1,NGEOM
       DO I=1,NCONV(IGEOM)
        X = VCONV(IGEOM,I)
        IF(X.LT.V1(1).OR.X.GT.V1(NP))THEN
         PRINT*,'Error in forwarderr.f : wavelength not covered'
         PRINT*,X,V1(1),V1(NP)
         STOP
        ENDIF
        CALL INTERP(V1,R1,NP,Y,X)
        RERR(IGEOM,I)=Y
       ENDDO
      ENDDO

      RETURN

      END
