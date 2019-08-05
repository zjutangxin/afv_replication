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
        implicit none
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
        implicit none

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

    subroutine PolynomialInt(xVec,nxR,yVec,nyR,zMx,zzMx,CoeVec)
    ! polynomial fit
        implicit none

        integer, intent(in) :: nxr,nyr
        real(dp), dimension(nxr), intent(in) :: xvec
        real(dp), dimension(nyr), intent(in) :: yvec
        real(dp), dimension(nxr,nyr), intent(in) :: zmx
        real(dp), dimension(nxr,nyr), intent(out) :: zzmx
        real(dp), dimension(7), intent(out) :: coevec

        integer :: indx,indy,indxy,info
        integer, parameter :: block_size = 32
        integer, parameter :: lwork = 7+7*block_size
        real(dp), dimension(nxr*nyr) :: DepVarVec
        real(dp), dimension(nxr*nyr,7) :: IndVarMx, IndVarMx_pre
        real(dp), dimension(lwork) :: work
        
        do indx = 1, nxR
            do indy = 1, nyR
                indxy=(indx-1)*nyR+indy
                IndVarMx(indxy,1) = 1.0_dp
                IndVarMx(indxy,2) = xVec(indx)
                IndVarMx(indxy,3) = yVec(indy)
                IndVarMx(indxy,4) = xVec(indx)*yVec(indy)
                IndVarMx(indxy,5) = xVec(indx)**2.0_dp
                IndVarMx(indxy,6) = yVec(indy)**2.0_dp
                IndVarMx(indxy,7) = xVec(indx)**2.0_dp*yVec(indy)**2.0_dp
                DepVarVec(indxy) = zMx(indx,indy)
            end do
        end do
        IndVarMx_pre = IndVarMx
        ! syntax
        ! DGELS(TRANS,M,N,NRHS,A,LDA,B,LDB,WORK,LWORK,INFO)
        call DGELS('N',nxr*nyr,7,1,IndVarMx,nxr*nyr,DepVarVec,nxr*nyr,work,lwork,info)
        coevec = DepVarVec
        ! CALL DRLSE(nxR*nyR,DepVarVec,6,IndVarMx(:,2:7),nxR*nyR,1,coeVec,SST,SSE)
        do indx = 1, nxR
            do indy = 1, nyR
                indxy = (indx-1)*nyR+indy
                zzMx(indx,indy) = dot_product(IndVarMx_pre(indxy,:),coevec)
            end do
        end do
    end subroutine PolynomialInt

    subroutine FindIndex(xx,xxVec,xxGrid,ii,pp)

        implicit none

        integer, intent(in) :: xxGrid
        integer, intent(out) :: ii
        real(dp), intent(in) :: xx
        real(dp), dimension(xxGrid), intent(in) :: xxVec
        real(dp), intent(out) :: pp

        if (xx <= xxVec(1)) then
            ii = 2
            pp = 1.0_dp
        elseif (xx >= xxVec(xxGrid)) then
            ii = xxGrid
            pp = 0.0_dp
        else
            ii = 1
            do while (xx > xxVec(ii))
                ii = ii+1
            end do
            pp = (xxVec(ii)-xx)/(xxVec(ii)-xxVec(ii-1))
        endif        
        
    end subroutine FindIndex

    function Afun(zz,pp)
    ! returns to capital
        use parameter, only: theta, zbar
        implicit none

        real(dp), intent(in):: zz, pp
        real(dp) :: afun

        afun = theta*zz/(zbar**(1-theta))+(zz-zbar)*pp/zbar
    end function Afun

    function phifun(pppr,bbpr)
    ! risky share
        use parameter, only: nz
        use global, only: zvec, pzvec
        implicit none

        real(dp), intent(in) :: pppr, bbpr
        real(dp) :: phifun
            
        integer :: indz
        real(dp) :: zzpr

        phifun = 0.0_dp
        do indz = 1, nz
            zzpr = zVec(indz)
            phifun = phifun + &
                ((Afun(zzpr,pppr)+pppr)/(Afun(zzpr,pppr)+pppr+bbpr)) &
                *pzVec(indz)
        end do
    end function phifun

    function rfun(bb,bbpr,pppr)
    ! bond return
        use parameter, only: zbar, bbeta
        implicit none

        real(dp), intent(in) :: bb,bbpr,pppr
        real(dp) :: rfun

        real(dp) :: pphi

        pphi = phifun(pppr,bbpr)
        rfun = (1-bbeta*pphi)*bbpr/(bbeta*(1-pphi)*(Afun(zbar,pppr)+bb))

    end function rfun

    subroutine piksr2(n,arr,brr)

        implicit none

        integer, intent(in) :: n
        real(dp), dimension(n), intent(inout) :: arr(n), brr(n)

        integer :: i,j
        real(dp) :: a,b

        DO  J=2,N
            A=ARR(J)
            B=BRR(J)
            DO I=J-1,1,-1
                IF (ARR(I) .LE. A) GO TO 10
                    ARR(I+1)=ARR(I)
                    BRR(I+1)=BRR(I)
            END DO
            I=0
            10    ARR(I+1)=A
            BRR(I+1)=B
        END DO        
    end subroutine piksr2

    function FindPerc(perc,a,p,n)
        ! Share of a by top perc in p
        implicit none

        integer, intent(in) :: n
        real(dp), dimension(n), intent(in) :: a,p
        real(dp), intent(in) :: perc
        real(dp) :: FindPerc

        real(dp), dimension(n) :: xcum,pcum
        integer :: indi
    
        pcum(1) = p(1)
        xcum(1) = a(1)*p(1)

        do indi = 2,n
            pcum(indi) = pcum(indi-1)+p(indi)
            xcum(indi) = xcum(indi-1)+a(indi)*p(indi)
        end do
        pcum = pcum/pcum(n)
        xcum = xcum/xcum(n)

        FindPerc = 1.0_dp - FFindInt(perc,pcum,xcum,n)

    end function FindPerc

    function FindGini(a,p,n)
        ! find Gini coefficient
        implicit none

        integer, intent(in) :: n
        real(dp), dimension(n), intent(in) :: a,p
        real(dp) :: FindGini

        real(dp), dimension(n) :: xcum,pcum,wwvec,hhvec
        integer :: indi

        pcum(1) = p(1)
        xcum(1) = a(1)*p(1)

        do indi = 2,n
            pcum(indi) = pcum(indi-1)+p(indi)
            xcum(indi) = xcum(indi-1)+a(indi)*p(indi)
        end do
        pcum = pcum/pcum(n)
        xcum = xcum/xcum(n)
        wwvec(1) = pcum(1)
        hhvec(1) = xcum(1)

        do indi = 2,n
            wwvec(indi) = pcum(indi) - pcum(indi-1)
            hhvec(indi) = xcum(indi) + xcum(indi-1)
        end do

        FindGini = 1.0_dp - dot_product(wwvec,hhvec)
    end function FindGini

end module routines