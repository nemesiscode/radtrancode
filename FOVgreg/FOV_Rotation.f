! FOV_Rotation is a Fortran 90/95 subroutine for computing the actual vertical FOV
!   given a rotation angle of the MCS focal planes.  
!
! Inputs:
!   Xfov(ih,ndet,narr) -- real    -- array of angular axes for horizontal FOV
!   Yfov(iv,ndet,narr) -- real    -- array of angular axes for vertical FOV
!   Hfov(ih,ndet,narr) -- real    -- array of horizontal FOV responses
!   Vfov(iv,ndet,narr) -- real    -- array of vertical FOV responses
!   HDet(ndet,narr)    -- real    -- array of horizontal detector positions
!   VDet(ndet,narr)    -- real    -- array of vertical detector positions
!   Lh(ndet,narr)      -- integer -- array of response lengths for horizontal FOV
!   Lv(ndet,narr)      -- integer -- array of response lengths for vertical FOV
!   ih                 -- integer -- maximum length of horizontal FOV response
!   iv                 -- integer -- maximum length of vertical FOV response
!   ndet               -- integer -- number of detectors
!   narr               -- integer -- number of arrays
!   angle              -- real    -- rotation angle, assumed to be in degrees
!   iy                 -- integer -- length of output arrays after interpolation
!                                    (note, iy = 2^N, where N is often 10)
!
! Outputs:
!   vgrid(iy)          -- real    -- regular angular grid for final vertical FOV
!   vout(iy,ndet,narr) -- real    -- final vertical FOV response
!
! Notes:
! * The FOV characterization is of that of the telescopes' focal plane as if one
!   is inside the instrument looking out.  The coordinates are (x,y) ==> (azimuth
!   angle,elevation angle).  
!   * The A arrays are ordered as follows: [A6 A5 A1 A2 A3 A4], with A6 having the 
!     most negative azimuth angle and A4 having the most positive azimuth angle.  
!   * The B arrays are ordered as follows: [B3 B2 B1], with B3 having the most 
!     negative azimuth angle and B1 having the most positive azimuth angle.
!   * The A detectors are ordered such that A_1 has the most negative elevation
!     angle and A_21 has the most positive elevation angle. 
!   * The B detectors are the opposite of the A detectors: i.e., B_1 has the most
!     positive elevation angle and B_21 has the most negative elevation angle.
!   * When mounted on MRO, MSC's elevation axis increases towards Mars.  That is,
!     detectors A_21 and B_1 are often viewing the surface of Mars, whereas the
!     other detectors, with lower elevation angles, are viewing the atmosphere.
! * To summarize the above, for Telescope A, the corner detectors are:
!     
!                          Negative Azimuth | Positive Azimuth
!                                           
!     Negative Elevation          A6_1      |      A4_1
!                            --------------------------------
!     Positive Elevation          A6_21     |      A4_21
!
!                               [ Surface of Mars here ]
!
! * The rotation angle is assumed to be input in degrees.  This is the rotation of
!   the fields of view about the boresight of the telescope (i.e., (0,0) in the 
!   (az,el) plane.  By convention, a POSITIVE ROTATION ANGLE IS CLOCKWISE.  To be 
!   clear about this, a positive rotation angle (for nominal fore-limb staring 
!   configuration) will lift A6_21 off of the surface of Mars and plunge A4_21
!   deeper down into the surface.  
!
! * For now, this subroutine will make use of FFT subroutines from Numerical 
!   Recipes that assume the array lengths involved are powers of two.  Hopefully
!   this can be upgraded later to be used with fftw ("the Fastest Fourier Transform
!   in the West") libraries so that it will NOT be restricted to powers of 2.  As
!   of now it appears the 2^10 (1024) is the target length for our arrays.  This 
!   allows for interpolation at resolution 2 mrad without chopping the edges of the
!   non-zero portions of the FOV curves.  
! * F90/95 allow for array inquiry functions that figure out array dimensions from 
!   what has been passed, but I believe it is faster to explicitly declare them...
!   and it is faster yet to have those declared sizes be parameters instead of 
!   variables.  So, we know that when using MCS data, ndet = 21 and narr = 9.  For
!   the first batch of FOV responses I have created, ih = 158 and iv = 190.  ih
!   and iv might change if the FOV curves are changed.  These values could be fed 
!   in through a module of parameters, or even a common block, if desired.
! * Each response curve has a potentially different length, so the arrays Lh and
!   Lv are used to refer only to the meaningful portions of the FOV arrays, e.g.,
!   Vfov(1:Lv(11,1),11,1) refers to the complete vertical FOV response for 
!   detector 11 on array A1.
! * The array numbering convention is NO LONGER:  [A4 A3 A2 A1 A5 A6 B1 B2 B3] = 
!                                                 [ 1  2  3  4  5  6  7  8  9]
!   But is instead now:  [A1 A2 A3 A4 A5 A6 B1 B2 B3] =
!                        [ 1  2  3  4  5  6  7  8  9]
! * There is hope that this subroutine can be modified to one day accept a 
!   specified grid on which to output the final vertical FOV responses, but for 
!   now we continue with assuming a grid of regular spacing of length iy = 2**N
! * This subroutine is written with expressions like "real" and 1.0e-1, which at 
!   first appear to be single precision -- if the user wants to run the program
!   in double precision, then this should be specified at compile time with a
!   "-r8" flag (or similar).
!
! This subroutine makes use of a routine we call interp.f90, which performs linear
!   interpolation.  interp.f90 is included below after the main subroutine in this 
!   file.
!
! This subroutine calls the subroutine convlv.f, which in turn calls realft.f and 
!   and twofft.f, both of which call four1.f.  convlv.f, realft.f, twofft.f, and 
!   four1.f are all from Numerical Recipes.  They have all been appended to the
!   end of this file. -- one note of experience from a user:  the Numerical 
!   Recipes routines make at least one "slight of hand" where a real array is 
!   passed but declared as complex upon receipt.  Everything should compile and 
!   run fine, but when a user tried to re-package this into a Fortran 90/95 
!   module, the compiler was not happy about the change of variable type.
!
!--------------------------------------------------------------------------------
! The approach of this subroutine is as follows:
! 1. We stretch the angular axes of the various detectors by the appropriate 
!    trigonometric scaling factors (cos & sin) to get the correct projections
!    onto the true vertical
! 2. We find the implied vertical angular translations for each detector owing
!    to the rotation
! 3. We interpolate the horizontal and vertical components on to regular grids
! 4. We perform a convolution of the two projections via the FFT method
!
!--------------------------------------------------------------------------------
! A further correction has become necessary to handle small but non-zero angles
!   of rotation (e.g., 0.1 degrees).  If the angle is small enough (or close 
!   enough to N.pi/2, N = integer), then the trigonometric projection of the
!   nearly perpendicular profile onto the vertical will be very small.  In 
!   practice, a focal plane rotation has two effects:  the detector positions
!   shift relative to their nominal (zero rotation angle) positions, and the 
!   shapes of the FOV profiles change via a convolution.  A small but non-zero
!   rotation angle can give a perceptible contribution for the former effect, 
!   while having no impact for the latter effect (because the perpendicular 
!   profile contributes essentially a delta function to the convolution).  If
!   one proceeds naively with small but non-zero rotation angles, then one can
!   incur various errors with the code, ranging from floating point errors to
!   just obtaining nonsensical FOV profiles.  Hence, we have implemented the 
!   following change:  If rotation angle is within 1.0 degrees of N.pi/2, then 
!   the convolution is no longer computed, however the shift of the detector
!   positions is still computed.  1.0 degrees was chosen because it is the 
!   angle where sin(theta) * the widest detector is approximately equal to the
!   interpolated grid spacing, 0.0002 radians.  4*(azimuthal half-width) of a
!   B3 detector (det #18) is used as the full width of the widest detector 
!   response signal (which is equal to 0.018576 radians), whereas 
!   4*(azimuthal half-width) of an A6 detector (det #20) is used as the full 
!   width of the narrowest detector response signal (which is 0.0026812 rad).
!   
!
! Thursday 20 July 2006, WGL
! Revised Friday 12 January 2007, WGL
! Re-revised Wednesday 7 February 2007, WGL

subroutine FOV_Rotation( Yfov, Xfov, Vfov, Hfov, VDet, HDet, Lv, Lh, iv, ih, &
                         ndet, narr, angle, iy, vgrid, vout )

  implicit none

  ! Subroutine arguments
  integer, intent(in)  :: ih, iv, ndet, narr, iy
  integer, intent(in), dimension(ndet,narr)    :: Lh, Lv
  real,    intent(in)  :: angle
  real,    intent(in), dimension(ndet,narr)    :: HDet, VDet
  real,    intent(in), dimension(ih,ndet,narr) :: Xfov, Hfov
  real,    intent(in), dimension(iv,ndet,narr) :: Yfov, Vfov
  real,    intent(out) :: vgrid(iy)
  real,    intent(out) :: vout(iy,ndet,narr)

  ! Local variables
  integer :: k, j, i, w, g, Ldat, ind(1)
  real    :: pi, pio2, rth, cosrth, sinrth, Vp(ndet,narr)  !, Rot(2,2)
  real    :: delt, adelt, thresh
  real, allocatable   :: xdat(:), ydat(:)
  real, dimension(iy) :: vcomp, hcomp, orig, resp
  real    :: ans(iy+iy), test1, testiy
  logical :: delta


  ! Deal with constants
  delt = 0.0002               ! 0.2 mrad resolution for the interpolation --
                              !   perhaps this should be made user-specified?
  pi = 4.0*atan(1.0)          ! pi may be defined in main program, so it could
                              !   possibly be passed or shared with this routine
  pio2 = pi/2.0
  adelt = delt*180.0/pi       ! delt in degrees
  thresh = 1.0                ! rotation angle of 1.0 degrees, used as a threshold 
                              !   for when to actually perform convolutions
  delta = .false.             ! logical switch to indicate whether one of the 
                              !   components is essentially a delta function as
                              !   far as the convolution is concerned

  ! Mathematically, a rotation matrix defines a positive rotation as counter-
  !   clockwise, but for MCS, a positive rotation angle is clockwise.  This is 
  !   not inconsistent since MCS's elevation axis (it's y-axis) is reverse from
  !   the traditional mathematical definition -- a positive rotation for a 
  !   traditional right-handed coordinate system is counter-clockwise, and a 
  !   positive rotation in MCS's "upside down" right-handed coordinate system
  !   is clockwise.  Therefore, the rotation angle sign remains unchanged.
  rth = angle*pi/180.0        ! angle in radians
  cosrth = cos(rth)
  sinrth = sin(rth)
  ! Rotation matrix = [cos(rth) -sin(rth); sin(rth) cos(rth)];
  ! We don't actually use the rotation matrix form for our conversions
  ! Rot(1,1) =  cosrth;  Rot(1,2) = -sinrth
  ! Rot(2,1) =  sinrth;  Rot(2,2) =  cosrth


  ! Construct the interpolated grid -- including the point 0 rad and centered 
  !   almost evenly about it
  vgrid(1) = delt*real( -iy/2 + 1 )
  do k = 2, iy
     vgrid(k) = vgrid(k-1) + delt
  end do

  array: do g = 1, 9
     detector: do w = 1, 21

        ! Find the vertical translation distance owing to the rotation --  
        !   this all depends on the defined locations of the detectors, which 
        !   is based on FWHM.  
        ! Full rotation should be found via linear algebra: x' = Rot*x,
        !   but we can skip this because we are only interested in the new
        !   vertical component:  vp = Rot(2,1)*HDet + Rot(2,2)*VDet
        ! Vp(w,g) = Rot(2,1)*HDet(w,g) + Rot(2,2)*VDet(w,g)  
        Vp(w,g) = sinrth*HDet(w,g) + cosrth*VDet(w,g) 

        ! If angle is 0, pi/2, pi, 3pi/2, or 2pi, then one component will not
        !   be used and its delta function contribution will break the 
        !   interpolation scheme.  Also, its delta function contribution will
        !   not change the convolution.

        ! Vertical FOV curve will not be used if angle = +/- pi/2 or +/- 3pi/2, 
        !   so test to see whether we are within the specified threshold of 
        !   either of those angles.
        ! Easier dealt with in degrees -- also, the routine needs to be able to
        !   handle negative angles, so we have made amends for that here.
        vert: if ( ( abs( abs(angle) - 90.0 ) < thresh ) .or. &
             ( abs( abs(angle) - 270.0 ) < thresh ) ) then

           vcomp = 0.0
           ! Need to set one element to unity -- Numerical Recipes requires
           !   "wrap around" format for one of the profiles to be convolved.
           ! Here we are sure that vcomp will be treated as response function
           !   instead of original function (in N.R. terms), so set last
           !   element to 1.
           vcomp(iy) = 1.0

           ! Set the logical switch delta to true
           delta = .true.

        ! If we are not in the special case, then we must perform linear 
        !   interpolation to find vcomp.  Add on zero end points so that the 
        !   interpolation is well-defined.  Scale and shift the angular axes.
        else

           ! Length of data curves to be interpolated (including added endpoints)
           Ldat = Lv(w,g)+2
           ! Allocate space for x and y data curves since Ldat is potentially
           !   different for each pair of (w,g)
           allocate( xdat(Ldat), ydat(Ldat) )

           ! Append angles farther separated than the grid onto which we want
           !   to interpolate -- these points have zero amplitude
           ! We need to make an monotonically increasing x-axis -- four cases
           !   to worry about: 
           !   1. cos*Yfov extends past vgrid limits and is increasing (th ~ 0)
           !   2. cos*Yfov within vgrid limits and is increasing (270 < th < 90)
           !   3. cos*Yfov within vgrid limits but is decreasing (90 < th < 270)
           !   4. cos*Yfov extends past vgrid limits but is decreasing (th ~ 180)
           test1 = cosrth*Yfov(1,w,g)
           testiy = cosrth*Yfov(Lv(w,g),w,g)
           
           ! Check to see which end has a larger angle (shouldn't be equal!)
           vmonotonic: if ( test1 < testiy ) then      ! monotically increasing

              ! Test if within vgrid lower limit 
              if ( test1 >= vgrid(1) ) then            ! is within
                 xdat(1) = vgrid(1)-delt
              else                                     ! not within
                 xdat(1) = test1-delt
              end if

              ! Test if within vgrid upper limit
              if ( testiy < vgrid(iy) ) then           ! is within
                 xdat(Ldat) = vgrid(iy)+delt
              else                                     ! not within
                 xdat(Ldat) = testiy+delt
              end if
              
              ! Fill in xdat and ydat
              xdat(2:Ldat-1) = (/ cosrth*( Yfov(:Lv(w,g),w,g)-VDet(w,g) ) + &
                                  Vp(w,g) /)
              ydat = (/ 0.0, Vfov(:Lv(w,g),w,g), 0.0 /)

           else                                        ! monotically decreasing
              
              ! Test if within vgrid lower limit 
              if ( testiy >= vgrid(1) ) then           ! is within
                 xdat(1) = vgrid(1)-delt
              else                                     ! not within
                 xdat(1) = testiy-delt
              end if

              ! Test if within vgrid upper limit
              if ( test1 < vgrid(iy) ) then            ! is within
                 xdat(Ldat) = vgrid(iy)+delt
              else                                     ! not within
                 xdat(Ldat) = test1+delt
              end if
              
              ! Fill in xdat and ydat -- consistently reverse angles and data to
              !   make monotonically increasing
              ydat(1) = 0.0;   ydat(Ldat) = 0.0
              do k = 1, Lv(w,g)
                 xdat(k+1) = cosrth*( Yfov(Lv(w,g)+1-k,w,g)-VDet(w,g) ) + Vp(w,g) 
                 ydat(k+1) = Vfov(Lv(w,g)+1-k,w,g)             
              end do

           end if vmonotonic

           ! Call our hand-written linear interpolation program
           call interpF( xdat, ydat, Ldat, vgrid, vcomp, iy )
           ! Deallocate the memory used by xdat and ydat, lest there be 
           !   memory leaks
           deallocate( xdat, ydat )

        end if vert

        ! Horizontal FOV curve will not be used if angle = 0, pi, or 2pi, so 
        !   test to see whether we are within the specified threshold of 
        !   either of those angles.
        ! Easier dealt with in degrees -- also, the routine needs to be able to
        !   handle negative angles, so we have made amends for that here.
        horiz: if ( ( abs( abs(angle) ) < thresh ) .or. &
                    ( abs( abs(angle) - 180.0 ) < thresh ) .or. &
                    ( abs( abs(angle) - 360.0 ) < thresh ) ) then

           hcomp = 0.0
           ! Need to set one element to unity -- Numerical Recipes requires
           !   "wrap around" format for one of the profiles to be convolved.
           ! Here we are sure that hcomp will be treated as response function
           !   instead of original function (in N.R. terms), so set last 
           !   element to 1.
           hcomp(iy) = 1.0

           ! Set the logical switch delta to true
           delta = .true.

        ! If we are not in the special case, then we must perform linear 
        !   interpolation to find hcomp.  Add on zero end points so that the 
        !   interpolation is well-defined.  Scale and shift the angular axes.
        else

           ! Length of data curves to be interpolated (including endpoints)
           Ldat = Lh(w,g)+2
           ! Allocate space for x and y data curves since Ldat is potentially
           !   different for each pair of (w,g)
           allocate( xdat(Ldat), ydat(Ldat) )

           ! Append angles farther separated than the grid onto which we want
           !   to interpolate -- these points have zero amplitude
           ! We need to make an monotonically increasing x-axis -- four cases
           !   to worry about: 
           !   1. sin*Xfov extends past vgrid limits and increasing (th = 0)
           !   2. sin*Xfov within vgrid limits and increasing (270 < th < 90)
           !   3. sin*Xfov within vgrid limits but decreasing (90 < th < 270)
           !   4. sin*Xfov extends past vgrid limits and decreasing (th = 180)
           test1 = sinrth*Xfov(1,w,g)
           testiy = sinrth*Xfov(Lh(w,g),w,g)
           
           ! Check to see which end has a larger angle (shouldn't be equal!)
           hmonotonic: if ( test1 < testiy ) then      ! monotically increasing

              ! Test if within vgrid lower limit 
              if ( test1 >= vgrid(1) ) then            ! is within
                 xdat(1) = vgrid(1)-delt
              else                                     ! not within
                 xdat(1) = test1-delt
              end if

              ! Test if within vgrid upper limit
              if ( testiy < vgrid(iy) ) then           ! is within
                 xdat(Ldat) = vgrid(iy)+delt
              else                                     ! not within
                 xdat(Ldat) = testiy+delt
              end if
              
              ! Fill in xdat and ydat
              xdat(2:Ldat-1) = (/ sinrth*( Xfov(:Lh(w,g),w,g)-HDet(w,g) ) + &
                                  Vp(w,g) /)
              ydat = (/ 0.0, Hfov(:Lh(w,g),w,g), 0.0 /)

           else                                        ! monotically decreasing
              
              ! Test if within vgrid lower limit 
              if ( testiy >= vgrid(1) ) then           ! is within
                 xdat(1) = vgrid(1)-delt
              else                                     ! not within
                 xdat(1) = testiy-delt
              end if

              ! Test if within vgrid upper limit
              if ( test1 < vgrid(iy) ) then            ! is within
                 xdat(Ldat) = vgrid(iy)+delt
              else                                     ! not within
                 xdat(Ldat) = test1+delt
              end if
              
              ! Fill in xdat and ydat -- consistently reverse angles and data to
              !   make monotonically increasing
              ydat(1) = 0.0;   ydat(Ldat) = 0.0
              do k = 1, Lh(w,g)
                 xdat(k+1) = sinrth*( Xfov(Lh(w,g)+1-k,w,g)-HDet(w,g) ) + Vp(w,g) 
                 ydat(k+1) = Hfov(Lh(w,g)+1-k,w,g)             
              end do

           end if hmonotonic

           ! Call our hand-written linear interpolation program
           call interpF( xdat, ydat, Ldat, vgrid, hcomp, iy )
           ! Deallocate the memory used by xdat and ydat, lest there be 
           !   memory leaks
           deallocate( xdat, ydat )

        end if horiz

        ! Now that the two components' contributions have been found, it is 
        !   time to compute their convolution.  To do this, we will use the 
        !   Numerical Recipes routines which operate via the FFT.  This is
        !   supposedly the most efficient method to accomplish this, but it
        !   requires some massaging of the data beforehand by the user to 
        !   make the routine happy.  In particular, one of the two curves is
        !   known as the "original function" and the other is known as the 
        !   "response function."  The response function is typically much
        !   shorter in length than the original function, so that is the main
        !   distinguishing feature.  The NR algorithm asks that the response
        !   function be stored in "wrap around" format, where the maximum value
        !   is in index 1 and the data points to the left of it are stored at 
        !   the upper end of the wrapped array:
        !
        !   _____---^--_____ becomes ^--_________---
        !
        ! Because the FFT is being used to compute the convolution, two 
        !   assumptions are being made (Numerical Recipes spells these out
        !   very clearly):  1. original and response functions have the same
        !   data length, and 2. both are periodic.  Our interpolation procedure
        !   solves the first potential problem -- many of the points are made
        !   to be zero.  The periodic problem is typically solved by zero 
        !   padding.  There are specific rules for handling this so that there
        !   is no periodic-wrap-around contamination, namely that the original
        !   function should be zero padded to the length equal to the number of
        !   response function points to the left of its maximum (i.e., the 
        !   points that have been wrapped around to the far end of the array).
        !   Since each detector's response is different, and further differences
        !   arise from each different rotation angle considered, we could (and 
        !   probably should) check to make sure each curve is sufficiently zero
        !   padded.  However, if we were to add more zeros, our data lengths 
        !   would have to increase to the next value of 2**N, and that's a 
        !   pretty big computational burden to alleviate potential leakage. 
        !   Also, since the interpolated functions already have a great deal of
        !   zero padding (even with rotation angle = 0 or 90), I'm not that 
        !   worried about it.  
        ! That said, there is probably potential to make things more efficient
        !   here by taking explicit advantage of all the zeros....

        ! Determine which component is the "response function" -- that is, which 
        !   component has more zero padding after the trigonometric rescaling?
        ! Probably the best way to do this is via searching for and counting the
        !   number of zeros in each component, but these operations could be 
        !   expensive when repeated over two components for each of 189 detectors.
        !   So, here we'll base the designation on angle.  If the rotation angle
        !   is within a certain number of degrees of pi/2 or 3pi/2, then we will
        !   designate the horizontal component as the original function and the 
        !   vertical component as response function.  The vertical component will
        !   be the original function in all other circumstances since they are 
        !   beginning with fewer zeros in their profiles.  From some testing
        !   with the profiles off-line, we choose the threshold angle to be:
        !   +/-26 degrees of 90 or 270.

        ! We will only do the convolution if our logical switch, delta, is still
        !   set to false.  Otherwise, there is no need to do the convolution
        !   because one function is acting as an effective delta function; 
        !   however, the detector position will have still moved due to the 
        !   rotation, so the final vertical profile will be determined by the
        !   interpolated component that isn't the delta function.
        delta1: if ( .NOT. delta ) then

           ! Horizontal component is original function & vertical component is 
           !   response
           wrap: if ( ( abs( angle - 90.0 ) < 26.0 ) .or. &
                ( abs( angle - 270.0 ) < 26.0 ) ) then
              
              orig = hcomp
           
              ! To conform to wrap-around format for the N.R. algorithm, find the 
              !   location of the maximum value of the response function
              ind = maxloc( vcomp )
              
              ! Now we need to shift all elements before this maximum to the end of
              !   the array -- we can use F90's circular shift function
              resp = cshift( vcomp, ind(1)-1 )

           ! Vertical component is original function & horiz. component is response
           else
           
              orig = vcomp
              ind = maxloc( hcomp )
              resp = cshift( hcomp, ind(1)-1 )

           end if wrap

           ! Now call Numerical Recipes' canned convolution routine
           call convlv(orig,iy,resp,iy,1,ans)

        ! If one component is a delta function, then the other one is the answer
        else

           ! Horizontal component is original function & vertical component is 
           !   response
           assign: if ( ( abs( angle - 90.0 ) < 26.0 ) .or. &
                ( abs( angle - 270.0 ) < 26.0 ) ) then
           
              ans = hcomp

           ! Vertical component is original function & horiz. component is response
           else

              ans = vcomp
              
           end if assign

        end if delta1

        ! The final vertical FOV is the first half of ans -- be sure to renormalize
        !   by max amplitude
        vout(:,w,g) = ans(:iy)/maxval(ans(:iy))

     end do detector
  end do array

end subroutine FOV_Rotation


!--------------------------------------------------------------------------------------
! This subroutine will perform linear interpolation from the input (xin,yin) onto a
!   specified output x-axis to yield (xout,yout).  xin and yin have dimension Nin, and
!   xout and yout have dimension Nout.  
! This subroutine assumes that the x axes are monotonically increasing, and it is 
!   crucial that they are -- effort is spent above in FOV_Rotation.f90 to ensure this.

subroutine interpF( xin, yin, Nin, xout, yout, Nout )

  implicit none

  integer, intent(in)  :: Nin, Nout
  real,    intent(in)  :: xin(Nin), yin(Nin), xout(Nout)
  real,    intent(out) :: yout(Nout)

  integer :: j, il, iu
  real    :: xl, xu, yl, yu, fac

  ! Beginning bounding index
  iu = 2

  do j = 1, Nout

     ! Check to see if we still have the right lower x index
     test: do 
        if ( xout(j) < xin(iu) ) then
           exit test
        else
           iu = iu + 1
        end if
     end do test
     il = iu - 1

     ! Designate bounding x and y values
     xl = xin(il);  xu = xin(iu)
     yl = yin(il);  yu = yin(iu)

     ! Linear interpolation factor
     fac = ( xout(j) - xl )/( xu - xl )

     yout(j) = yl + fac*( yu - yl )

  end do

end subroutine interpF


!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! Below are algorithms from the Numerical Recipes collection
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

! This file holds the pertinent Numerical Recipes files for doing convolution
!   via the FFT.  Minimal changes have been made to make them F90 compatible (i.e.,
!   comments changed from "c" to "!" and line continuation from 6th column to "&")
!
! Contents:
!   convlv.f
!   twofft.f
!   realft.f
!   four1.f

!-----------------------------------------------------------------------
      SUBROUTINE convlv(data,n,respns,m,isign,ans)
      INTEGER isign,m,n,NMAX
      REAL data(n),respns(n)
      COMPLEX ans(n)
      PARAMETER (NMAX=4096)
!U    USES realft,twofft
      INTEGER i,no2
      COMPLEX fft(NMAX)
      do 11 i=1,(m-1)/2
        respns(n+1-i)=respns(m+1-i)
11    continue
      do 12 i=(m+3)/2,n-(m-1)/2
        respns(i)=0.0
12    continue
      call twofft(data,respns,fft,ans,n)
      no2=n/2
      do 13 i=1,no2+1
        if (isign.eq.1) then
          ans(i)=fft(i)*ans(i)/no2
        else if (isign.eq.-1) then
          if (abs(ans(i)).eq.0.0) pause &
               'deconvolving at response zero in convlv'
          ans(i)=fft(i)/ans(i)/no2
        else
          pause 'no meaning for isign in convlv'
        endif
13    continue
      ans(1)=cmplx(real(ans(1)),real(ans(no2+1)))
      call realft(ans,n,-1)
      return
      END

!-----------------------------------------------------------------------
      SUBROUTINE twofft(data1,data2,fft1,fft2,n)
      INTEGER n
      REAL data1(n),data2(n)
      COMPLEX fft1(n),fft2(n)
!U    USES four1
      INTEGER j,n2
      COMPLEX h1,h2,c1,c2
      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5)
      do 11 j=1,n
        fft1(j)=cmplx(data1(j),data2(j))
11    continue
      call four1(fft1,n,1)
      fft2(1)=cmplx(aimag(fft1(1)),0.0)
      fft1(1)=cmplx(real(fft1(1)),0.0)
      n2=n+2
      do 12 j=2,n/2+1
        h1=c1*(fft1(j)+conjg(fft1(n2-j)))
        h2=c2*(fft1(j)-conjg(fft1(n2-j)))
        fft1(j)=h1
        fft1(n2-j)=conjg(h1)
        fft2(j)=h2
        fft2(n2-j)=conjg(h2)
12    continue
      return
      END

!-----------------------------------------------------------------------
      SUBROUTINE realft(data,n,isign)
      INTEGER isign,n
      REAL data(n)
!U    USES four1
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      END

!-----------------------------------------------------------------------
      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
