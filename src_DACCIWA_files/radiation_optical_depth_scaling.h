! radiation_optical_depth_scaling.h - Cloud optical-depth scaling for Tripleclouds 
!
! Copyright (C) 2016 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! This file is intended to be included inside a module to ensure that
! this simple routine may be inlined

!---------------------------------------------------------------------
pure subroutine optical_depth_scaling(nreg, frac_std, od_scaling)

  use parkind1, only : jprb

  ! Number of regions
  integer, intent(in)     :: nreg

  ! Fractional standard deviation of in-cloud water content
  real(jprb), intent(in)  :: frac_std

  ! Optical depth scaling for the cloudy regions
  real(jprb), intent(out) :: od_scaling(2:nreg)

  if (nreg == 2) then
    ! Only one clear-sky and one cloudy region: cloudy region is
    ! homogeneous
    od_scaling(2) = 1.0_jprb
  else
    ! Two cloudy regions with optical depth scaled by 1-x and
    ! 1+x.
    ! Simple version which fails when fractional_std >= 1:
    !od_scaling(2) = 1.0_jprb-cloud%fractional_std(jcol,jlev)
    ! According to Shonk and Hogan (2008), 1-x should correspond to the
    ! 16th percentile. If we treat this as a lognormal such that the
    ! equivalent Normal has a mean mu and standard deviation sigma, then
    ! the 16th percentile of the lognormal is very close to
    ! exp(mu-sigma).
    od_scaling(2) &
         &  = exp(-sqrt(log(frac_std**2+1))) / sqrt(frac_std**2+1)
    ! Ensure mean optical depth is conserved
    od_scaling(3) = 2.0_jprb-od_scaling(2)
  end if

end subroutine optical_depth_scaling
