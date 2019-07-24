! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Main subroutine implementing Property 1
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
! =============================================================================

subroutine SolveSystem(t,b1,b2,b1pr,b2pr,retval)

    use parameter
    use global
    use routines
    implicit none

    integer, intent(in) :: t
    real(dp), intent(in) :: b1,b2,b1pr,b2pr
    real(dp), dimension(14), intent(out) :: retval
    real(dp), parameter :: x_min = 1e-9

    real(dp) :: be1,be1pr,p1,p1pr,phi1,a1,r1,c1_e,EU1_e,c1_w,U1_w
    real(dp) :: be2,be2pr,p2,p2pr,phi2,a2,r2,c2_e,EU2_e,c2_w,U2_w
    real(dp) :: v1wpr, v2wpr, v1epr, v2epr
    real(dp) :: alfa_t, eta_t
    integer :: indz

    if (regime == 0) then
        be1   = wgt1*b1
        be2   = wgt2*b2
        be1pr = wgt1*b1pr
        be2pr = wgt2*b2pr
    else
        be1   = sbar*wgt1*(b1+b2)
        be2   = (1-sbar)*wgt2*(b1+b2)
        be1pr = sbar*wgt1*(b1pr+b2pr)
        be2pr = (1-sbar)*wgt2*(b1pr+b2pr)
    endif

    alfa_t = (1.0_dp-bbeta**(nT-t+1))/(1.0_dp-bbeta)
    eta_t  = (bbeta-bbeta**(nT-t+1))/(1.0_dp-bbeta**(nT-t+1))

    p1pr = AppFun2D(bVec,Nb,bVec,Nb,PriMx1,b1pr,b2pr)
    p2pr = AppFun2D(bVec,Nb,bVec,Nb,PriMx2,b1pr,b2pr)
    phi1 = phifun(p1pr,be1pr)
    phi2 = phifun(p2pr,be2pr)
	p1   = eta_t*phi1*(Afun(zbar,p1pr)+be1)/(1.0-eta_t*phi1)
    p2   = eta_t*phi2*(Afun(zbar,p2pr)+be2)/(1.0-eta_t*phi2)   
    
    if (t == nT) then
        R1    = 1.0_dp
        R2    = 1.0_dp
        c1_e  = Afun(zbar,p1) + be1
        c2_e  = Afun(zbar,p2) + be2
        EU1_e = log(1.0_dp-eta_t)
        EU2_e = log(1.0_dp-eta_t)
        do indz=1, Nz
            a1 = max(Afun(zVec(indz),p1)+p1+be1,x_min)
            a2 = max(Afun(zVec(indz),p2)+p2+be2,x_min)
            EU1_e = EU1_e + log(a1)*pzVec(indz)
            EU2_e = EU2_e + log(a2)*pzVec(indz)
        end do
    else
        R1    = (1.0_dp-eta_t*phi1)*be1pr/(eta_t*(1-phi1)*(Afun(zbar,p1)+be1))
        R2    = (1.0_dp-eta_t*phi2)*be2pr/(eta_t*(1-phi2)*(Afun(zbar,p2)+be2))
        c1_e  = ((1.0_dp-eta_t)/eta_t)*(p1+be1pr/R1)
        c2_e  = ((1.0_dp-eta_t)/eta_t)*(p2+be2pr/R2)
        EU1_e = log(1-eta_t)+(alfa_t-1)*log(eta_t*phi1/p1)
        EU2_e = log(1-eta_t)+(alfa_t-1)*log(eta_t*phi2/p2)
        do indz=1, Nz
            a1 = max(Afun(zVec(indz),p1)+p1+be1,x_min)
            a2 = max(Afun(zVec(indz),p2)+p2+be2,x_min)
            EU1_e = EU1_e+alfa_t*log(a1)*pzVec(indz)
            EU2_e = EU2_e+alfa_t*log(a2)*pzVec(indz)
        end do
    endif

    c1_w = wbar+wgt1*(b1pr/R1-b1)
    c2_w = wbar+wgt2*(b2pr/R2-b2)
    U1_w = log(max(c1_w,x_min))
    U2_w = log(max(c2_w,x_min))

    V1wpr = AppFun2D(bVec,Nb,bVec,Nb,v1wMx,b1pr,b2pr)
    V2wpr = AppFun2D(bVec,Nb,bVec,Nb,v2wMx,b1pr,b2pr)
    V1epr = AppFun2D(bVec,Nb,bVec,Nb,v1eMx,b1pr,b2pr)
    V2epr = AppFun2D(bVec,Nb,bVec,Nb,v2eMx,b1pr,b2pr)

    !----CHECK IF PRICES ARE POSITIVE----!

    if (p1 < 0) then
        EU1_e = -10000.0_dp
        U1_w  = -10000.0_dp
    endif
    if (p2 < 0) then
        EU2_e = -10000.0_dp
        U2_w  = -10000.0_dp
    endif    

    ! send the equilibrium values back to main
    retval(1)  = c1_e
    retval(2)  = c1_w
	retval(3)  = p1
	retval(4)  = R1
    retval(5)  = U1_w + delta*V1wpr
    retval(6)  = EU1_e + bbeta*V1epr
	retval(7)  = c2_e
	retval(8)  = c2_w
	retval(9)  = p2
	retval(10) = R2
    retval(11) = U2_w+delta*V2wpr
    retval(12) = EU2_e+bbeta*V2epr
    retval(13) = p1pr
    retval(14) = p2pr    

    contains
    
        function Afun(zz,pp)
        ! returns to capital
            real(dp), intent(in):: zz, pp
            real(dp) :: afun

            afun = theta*zz/(zbar**(1-theta))+(zz-zbar)*pp/zbar
        end function Afun

        function phifun(pppr,bbpr)
        ! risky share
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

end subroutine SolveSystem