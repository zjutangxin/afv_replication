! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Module with shared variables.
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

module global

    use parameter
    implicit none

    real(dp), dimension(nb) :: bvec
    real(dp), dimension(nz) :: zvec, pzvec
    real(dp), dimension(nb,nb) :: v1wMx, v2wMx, v1eMx, v2eMx, priMx1, priMx2
    real(dp), dimension(nb,nb) :: DebPol1,DebPol2
    real(dp), dimension(maxgrid) :: DebChoiceVec,respfun1,respfun2
    real(dp), dimension(mgrid) :: adist

end module global