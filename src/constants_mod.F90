!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief Define various constants for the application.
!>
!> @details Various computational, physical and geometrical constants are
!>          defined in this module. Their values are also set here.
module constants_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64

  implicit none

  ! Define default application-defined kinds for all intrinsic data types

  !> @name Set up default kinds for real and double-precision variables.
  !> @{
  real,             private :: r_val  !< A native real used to compute kind of native real.
  double precision, private :: dp_val !< A native double-precision used to compute kind of native dp.

  integer, parameter :: r_def     = real64 !< Default real kind for application.
  integer, parameter :: r_single  = real32 !< Default single precision real kind for application.
  integer, parameter :: r_double  = real64 !< Default double precision real kind for application.

  integer, parameter :: r_native  = kind(r_val)  !< Native kind for real.
  integer, parameter :: dp_native = kind(dp_val) !< Native kind for double precision.

  ! Define kinds specifically for IO
  integer, parameter :: dp_xios = kind(dp_val) !< XIOS kind for double precision fields

  integer, parameter :: r_um = real64 !< Default real kind used by the UM.

  !> @}

  !> @name Complex
  !> @{
  !> @}

  !> @name Set up default kinds for integers.
  !> @{
  integer, private   :: i_val     !< A native integer used to compute kind of native integer.

  integer, parameter :: i_def        = int32       !< Default integer kind for application.
  integer, parameter :: i_native     = kind(i_val) !< Native kind for integer.
  integer, parameter :: i_long       = int64       !< Explicit long integer.
  integer, parameter :: i_halo_index = int64       !< Integer kind for the index used in halo swapping
  integer, parameter :: i_um         = int32       !< Default integer kind used by the UM.
  !> @}

  !> @name Set up default kinds for logicals.
  !> @{
  logical, private   :: l_val     !< A native logical used to compute kind of native logical.

  integer, parameter :: l_def     = kind(l_val) !< Default logical kind for application.
  integer, parameter :: l_native  = kind(l_val) !< Native kind for logical.
  !> @}

  !> @name Set up default kinds for character variables.
  !> @{
  character, private :: c_val     !< A native character used to compute kind of native character.

  integer, parameter :: c_def     = kind(c_val) !< Default character kind for application.
  integer, parameter :: c_native  = kind(c_val) !< Native kind for character.
  !> @}

  !> @name Set up default lengths for string variables.
  !> @{
  integer, parameter :: str_short        = 16  !< Length of "short" strings.
  integer, parameter :: str_def          = 128 !< Default string length for normal strings.
  integer, parameter :: str_long         = 255 !< Default length of long string.
  integer, parameter :: str_max_filename = 512 !< Default maximum length of a file-name.
  !> @}

  !> @name Platform constants
  !> @{
  real(kind=r_def), parameter :: LARGE_REAL_POSITIVE = huge(0.0_r_def) !< The largest
  !<                            positive number of kind r_def that is not an infinity.
  real(kind=r_def), parameter :: LARGE_REAL_NEGATIVE = -LARGE_REAL_POSITIVE !< The largest
  !<                            negative number of kind r_def that is not an infinity.
  integer, parameter :: cache_block = 256 !< Size of a cache block, for padding
  !<                           arrays to ensure access to different cache lines

  !> @}

  !> @name Numerical constants
  !> @{
  real(kind=r_def), parameter :: EPS = 3.0e-15_r_def !< Relative precision: if (abs(x-y) > EPS) then assume x==y.
  !> @}

  !> @name Mathematical constants
  !> @{
  real(kind=r_def), parameter :: PI  = 4.0_r_def*atan(1.0_r_def) !< Value of pi.
  !> @}

  !> @name Conversion factors
  !> @{
  real(r_def), parameter :: degrees_to_radians = PI / 180.0_r_def
  real(r_def), parameter :: radians_to_degrees = 180.0_r_def / PI

  !> @}
  ! Missing data indicators
  real    (r_def), parameter :: RMDI = -huge(0.0_r_def)        !< Missing data indicator value for real numbers
  integer (i_def), parameter :: IMDI = -huge(0_i_def)          !< Missing data indicator value for integer numbers
  !> @}

end module constants_mod

