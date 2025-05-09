! This file contains the two subroutines necessary to read-in the FOV calibration
!   data.  At the moment, these routines assume that the FOV calibration data is
!   stored in a file called "Master_FOV.dat", which resides in the local directory
!   whence the program is executed.  
!
! The reading-in process is split into two steps:
!   1. First we read the data file and examine it to come up with array bound
!      information that may change from time to time, if the FOV data is ever
!      updated.  We are after the maximum data lengths for the FOV response
!      curves.  This comes from finding the maximum values of array lengths.
!   2. After we have the appropriate array lengths, then we can allocate array
!      space in which to store our angular axes and response curves.
!
! These routines assume that Master_FOV.dat has 189 columns, one for each MCS 
!   detector.
!
! Master.dat has a number of rows stored in the following order of sections:
!   1  Array Label (1 = A, 2 = B)
!   2  Channel Label (1-6, 1-3)
!   3  Detector Number (1-21)
!   4  Vertical detector position
!   5  Length of actual data stored in real arrays
!   6  Angular elevation axis of vertical FOV data
!   7  Vertical FOV response data
!   8  Horizontal detector position
!   9  Length of actual data stored in real arrays
!   10 Angular azimuth axis of horizontal FOV data
!   11 Horizontal FOV response data
!
! Sections 5 and 9 are the sections that give us array dimension size.  
! Sections 6 & 7 will have a length equal the maximum value of section 5, and 
!   any detectors which have fewer data points specified than this max value 
!   will have been zero padded at their end.  Similarly, Sections 10 & 11 will
!   have a length equal to the maximum value of section 9.
! Master_FOV.dat was generated by a Matlab script "ToWriteFOV.m" that reads in
!   the native Matlab saved variables in FinalVFOV.mat and FinalHFOV.mat.
!
!---------------------------------------------------------------------------------

! FOV_ReadIn_Sizes should be called first to retrieve the dimensions iv & ih.  
!   This takes in the dimension ndet and narr, which could potentially be hard-
!   coded into the routine or passed via a module or common block, and reads the
!   file "Master_FOV.dat" which holds all the FOV calibration data.
! The data in Master_FOV.dat were written in double precision, so this subroutine
!   forces one to assume this.

subroutine FOV_ReadIn_Sizes( ndet, narr, iv, ih )

  implicit none

  integer, intent(in)  :: ndet, narr
  integer, intent(out) :: iv, ih

  integer :: k 
  real*8  :: rdum(ndet*narr)
  character*100 aname

  ! Open the file containing the FOV calibration data -- ascii text, double
  !   precision storage

  ! Look for FOV file in raddata/
  aname='MCS_Master_FOV.dat'
  call datarchive(aname)  
  open(unit=10,file=aname,status='old')

  ! Read one row at a time
  !   The fifth row accesses Section 5 of the file, and this is what we want
  do k = 1, 5
     read(10,*) rdum
  end do

  ! Find the maximum value of this section, and convert it to an integer
  iv = nint( maxval( rdum ) )

  ! Using the array length just learned, read past the rows contained in 
  !   Sections 6 & 7 & 8
  do k = 1, iv+iv+1
     read(10,*) rdum
  end do

  ! Now read Section 9, containing the lengths of horizontal response data
  read(10,*) rdum

  ! Find the maximum value of this section, and convert it to an integer
  ih = nint( maxval( rdum ) )

  close(10)

end subroutine FOV_ReadIn_Sizes

!---------------------------------------------------------------------------------
! FOV_ReadIn should be called after memory has been allocated for the FOV arrays,
!   namely, Xfov, Yfov, Hfov, and Vfov.  These arrays depend on knowing ih & iv, 
!   which are found by calling FOV_ReadIn_Sizes.
!
! This subroutine reopens the FOV calibration data file, "Master_FOV.dat", and now
!   stores its contents in different arrays.  
!
! Inputs:
!   ndet               -- integer -- number of detectors
!   narr               -- integer -- number of arrays
!   iv                 -- integer -- maximum length of vertical FOV response
!   ih                 -- integer -- maximum length of horizontal FOV response
!
! Outputs:
!   Yfov(iv,ndet,narr) -- real    -- array of angular axes for vertical FOV
!   Xfov(ih,ndet,narr) -- real    -- array of angular axes for horizontal FOV
!   Vfov(iv,ndet,narr) -- real    -- array of vertical FOV responses
!   Hfov(ih,ndet,narr) -- real    -- array of horizontal FOV responses
!   VDet(ndet,narr)    -- real    -- array of vertical detector positions
!   HDet(ndet,narr)    -- real    -- array of horizontal detector positions
!   Lv(ndet,narr)      -- integer -- array of response lengths for vertical FOV
!   Lh(ndet,narr)      -- integer -- array of response lengths for horizontal FOV

subroutine FOV_ReadIn( ndet, narr, iv, ih, Yfov, Xfov, Vfov, Hfov, VDet, HDet, Lv, Lh)

  implicit none

  integer, intent(in)  :: ndet, narr, iv, ih
  integer, intent(out), dimension(ndet,narr)    :: Lv, Lh
  real,    intent(out), dimension(ndet,narr)    :: VDet, HDet
  real,    intent(out), dimension(iv,ndet,narr) :: Yfov, Vfov
  real,    intent(out), dimension(ih,ndet,narr) :: Xfov, Hfov

  integer :: k
  real*8  :: rdum(ndet*narr)

  ! Open the file containing the FOV calibration data -- ascii text, double
  !   precision storage
  open(unit=10,file='/home/jupiter/plan/teanby/mcs/FOVgreg/Master_FOV.dat',status='old')

  ! Read one row at a time
  ! Scroll through detector header label (Sections 1 - 3)
  do k = 1, 4
     read(10,*) rdum
  end do

  ! Read Section Vertical detector position (Section 4)
  VDet = reshape( rdum, (/ ndet, narr /) )

  ! Read Lengths of vertical response curves (Section 5), and convert to integers
  read(10,*) rdum
  Lv = nint( reshape( rdum, (/ ndet, narr /) ) )

  ! Read Vertical (elevation) angular axes (Section 6)
  do k = 1, iv
     read(10,*) rdum
     Yfov(k,:,:) = reshape( rdum, (/ ndet, narr /) )
  end do

  ! Read Vertical response curves (Section 7)
  do k = 1, iv
     read(10,*) rdum
     Vfov(k,:,:) = reshape( rdum, (/ ndet, narr /) )
  end do

  ! Read Horizontal detector position (Section 8)
  read(10,*) rdum
  HDet = reshape( rdum, (/ ndet, narr /) )

  ! Read Lengths of horizontal response curves (Section 9), and convert to integers 
  read(10,*) rdum
  Lh = nint( reshape( rdum, (/ ndet, narr /) ) )

  ! Read Horizontal (azimuth) angular axes (Section 10)
  do k = 1, ih
     read(10,*) rdum
     Xfov(k,:,:) = reshape( rdum, (/ ndet, narr /) )
  end do

  ! Read Horizontal response curves (Section 11)
  do k = 1, ih
     read(10,*) rdum
     Hfov(k,:,:) = reshape( rdum, (/ ndet, narr /) )
  end do
  
  close(10)

end subroutine FOV_ReadIn
