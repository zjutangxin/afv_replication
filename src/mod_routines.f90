! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Module with auxilliary routines.
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

module routines
    implicit none
    integer, parameter, private :: dp = kind(1.0d0)
    contains

    function FFindInt(xx,xxVec,yyVec,nn)
    ! one-dimension interpolation
        integer, intent(in) :: nn

        real(dp), intent(in) :: xx
        real(dp), dimension(nn) :: xxvec, yyvec
        real(dp) :: FFindInt

        integer :: ind
        real(dp) :: slope

        if (xx <= xxVec(1)) then
            slope    = (yyVec(2)-yyVec(1))/(xxVec(2)-xxVec(1))
            FFindInt = yyVec(1)+slope*(xx-xxVec(1))
        elseif (xx >= xxVec(nn)) then
            slope    = (yyVec(nn)-yyVec(nn-1))/(xxVec(nn)-xxVec(nn-1))
            FFindInt = yyVec(nn)+slope*(xx-xxVec(nn))
        else
            ind = 1
            do while (xx > xxVec(ind))
                ind = ind+1
            end do
            slope    = (yyVec(ind)-yyVec(ind-1))/(xxVec(ind)-xxVec(ind-1))
            FFindInt = yyVec(ind-1)+slope*(xx-xxVec(ind-1))
        endif
    end function FFindInt    

    function AppFun2D(xVec,nxR,yVec,nyR,zMx,x,y)
    ! two-dimension interpolation
        integer, intent(in) :: nxr, nyr
        real(dp), dimension(nxr), intent(in) :: xvec
        real(dp), dimension(nyr), intent(in) :: yvec
        real(dp), dimension(nxr,nyr), intent(in) :: zmx
        real(dp), intent(in) :: x, y
        real(dp) :: appfun2d
        
        integer :: nx,ny
        real(dp) :: xstep,ystep,px,py,z1,z2
      
        !------FIND LOCATION OF POINT x------!
      
        xstep = (xVec(nxR)-xVec(1))/(nxR-1)
        if (x <= xVec(1)) then
            nx = 1
            px = 1.0_dp
        elseif (x >= xVec(nxR)) then
            nx = nxR-1
            px = 0.0_dp
        else
            nx = 1 + floor((x-xVec(1))/xstep)
            px = 1.0_dp - (x-xVec(nx))/xstep
        endif
      
        !------FIND LOCATION OF POINT y------!
      
        ystep = (yVec(nyR)-yVec(1))/(nyR-1)
        if (y<=yVec(1)) then 
            ny = 1
            py = 1.0_dp
        elseif (y >= yVec(nyR)) then
            ny = nyR-1
            py = 0.0_dp
        else
            ny = 1 + floor((y-yVec(1))/ystep)
            py = 1.0_dp - (y-yVec(ny))/ystep
        endif
      
        !------COMPUTE z------!
      
        z1 = zMx(nx,ny)*px + zMx(nx+1,ny)*(1.0_dp-px)
        z2 = zMx(nx,ny+1)*px + zMx(nx+1,ny+1)*(1.0_dp-px)
        appfun2d = z1*py + z2*(1.0_dp-py)
      
    end function AppFun2D    

end module routines