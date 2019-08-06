! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Module with parameters shared across files.
! Author:
!     Xin Tang @ IMF, Summer 2019
!  
! Record of Revisions:
!         Date:                 Description of Changes
!     ===========        =================================
!      07/23/2019:          Original Code: No paralleling
!
! Compiling Environment:
!   GNU gfortran on Ubuntu 16.04
!
! Library Used:
!   - MINPACK: by source code
!   - LAPACK, ATLAS, and BLAS: by binary library
! 
! Shared by:
!   - debt_main.f90
! =============================================================================

module parameter

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    real(dp), parameter :: theta = 0.2_dp
    real(dp), parameter :: zbar = 1.0_dp, zmin = 0.86_dp
    real(dp), parameter :: zmax = 2.0_dp*zbar - zmin        
    real(dp), parameter :: wbar = (1-theta)*zbar**theta
    real(dp), parameter :: wgt1 = 0.796_dp, wgt2 = 0.796_dp
    real(dp), parameter :: sbar = 0.5_dp
    real(dp), parameter :: surv_rate = 0.975_dp
    real(dp), parameter :: rbar = 0.03_dp
    real(dp), parameter :: bbeta = surv_rate/(1+rbar) 
    real(dp), parameter :: delta = bbeta

    integer, parameter :: nb = 20, maxGrid = 200
    integer, parameter :: nt = 20
    integer, parameter :: regime = 1

    integer, parameter :: nz = 21
    real(dp), parameter :: bmin = 0.0_dp, bmax = 0.8_dp
    real(dp), parameter :: bbar = (bmin+bmax)/2.0_dp

    integer, parameter :: mgrid = 5000
    real(dp), parameter :: amin = 1e-7, amax = 1000.0_dp

    integer, parameter :: root = 0
    integer :: myrank

end module parameter