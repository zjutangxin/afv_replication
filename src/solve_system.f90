! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - A collection of major subroutines:
!       1. SolveSystem: implements Property 1 in the Appendix
!       2. NashSolution: Solve the Nash Equilibrium
!       3. Find_ss: Find the implied steady state bond
!       4. ss_distribution: find the stationary distribution of entrepreneurs
!       5. t_step_distribution: simulates the distribution forward
!
! Author:
!     Xin Tang @ IMF, Summer 2019
!     Based on Vincenzo Quadrini's original code
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

end subroutine SolveSystem

subroutine NashSolution(nn,xx,fvec,iflag)

    use parameter, only: maxgrid,dp
    use global, only: respfun1,respfun2,DebChoiceVec
    use routines
    implicit none

    integer, intent(in) :: nn
    integer, intent(out) :: iflag
    real(dp), dimension(nn), intent(in) :: xx
    real(dp), dimension(nn), intent(out) :: fvec

    real(dp) :: bb1pr, bb2pr
    
    bb2pr = FFindInt(XX(1),DebChoiceVec,RespFun2,MaxGrid)
    bb1pr = FFindInt(XX(2),DebChoiceVec,RespFun1,MaxGrid)
    fvec(1) = XX(1)-bb1pr
    fvec(2) = XX(2)-bb2pr

    iflag = 0

end subroutine NashSolution

subroutine find_ss(nn,xx,fvec,iflag)

    use parameter, only: dp, nb
    use global, only: bvec, DebPol1, DebPol2
    use routines
    implicit none

    integer, intent(in) :: nn
    integer, intent(out) :: iflag
    real(dp), dimension(nn), intent(in) :: xx
    real(dp), dimension(nn), intent(out) :: fvec
    real(dp) :: bb1,bb2,bb1pr,bb2pr

    bb1 = XX(1)
    bb2 = XX(2)
    bb1pr = AppFun2D(bVec,Nb,bVec,Nb,DebPol1,bb1,bb2)
    bb2pr = AppFun2D(bVec,Nb,bVec,Nb,DebPol2,bb1,bb2)
    fvec(1) = bb1pr-bb1
    fvec(2) = bb2pr-bb2

    iflag = 0

end subroutine find_ss

subroutine ss_distribution(b,p,Mea,stats)

    use parameter
    use global, only: zvec,pzvec,adist
    use routines

    implicit none
    real(dp), intent(in) :: b,p
    real(dp), dimension(mgrid), intent(out) :: Mea
    real(dp), dimension(3), intent(out) :: stats

    real(dp) :: w_size,be,phi,R,a,apr,wInc,perc
    real(dp), dimension(mgrid) :: MMea
    real(dp), dimension(mgrid*nz) :: inc_vec,mea_vec
    real(dp), dimension(mgrid,nz) :: prodist,ydist
    integer, dimension(mgrid,nz) :: inddist
    integer :: indi, indj, meaiter
    real(dp) :: MeaErr
    integer :: index, index_bar
    real(dp) :: prob, prob_bar

    w_size = wgt1

    be = w_size*b
    phi = phifun(p,be)
    R = Rfun(be,be,p)

    mea = 1.0_dp/real(mGrid)

    do indi = 1,mgrid
        a = adist(indi)
        do indj = 1,nz
            apr = bbeta*((Afun(zVec(indj),p)+p)*phi/p+R*(1-phi))*a
            call FindIndex(apr,aDist,mgrid,index,prob)
            inddist(indi,indj) = index
            prodist(indi,indj) = prob
        end do
    end do

    ! new entrants
    apr = afun(zbar,p)+p+be
    call findindex(apr,adist,mgrid,index_bar,prob_bar)

    ! iterating on wealth distribution
    MeaErr = 1000.0_dp
    MeaIter = 1

    do while ((MeaErr>1e-15) .and. (MeaIter<15000))
        MMea = 0.0_dp
        do indi = 1, mgrid
            do indj = 1, nz
                index = inddist(indi,indj)
                prob = prodist(indi,indj)
                MMea(index-1) = MMea(index-1)+&
                    surv_rate*Mea(indi)*prob*pzvec(indj)
                MMea(index) = MMea(index)+&
                    surv_rate*Mea(indi)*(1-prob)*pzvec(indj)
            end do
        end do
        MMea(index_bar-1) = MMea(index_bar-1)+(1-surv_rate)*prob_bar
        MMea(index_bar) = MMea(index_bar)+(1-surv_rate)*(1-prob_bar)
        MeaErr = sum(abs(MMea-Mea))

        if (mod(Meaiter,250)==0 .and. myrank == root) then
            write (*,'(A15,I5,ES14.6)') 'Dist Iter = ', Meaiter, Meaerr
        endif

        Mea = MMea
        MeaIter = MeaIter + 1
    end do
    Mea = Mea/sum(Mea)

    ! construct income distribution
    do indi = 1,mgrid
        a = adist(indi)
        ! ybardist(indi) = bbeta*(Afun(zbar,p)*phi/p+(R-1)*(1-phi))*a
        do indj = 1,nz
            ydist(indi,indj) = bbeta*(Afun(zVec(indj),p)*phi/p + &
                (R-1)*(1-phi))*a
        end do
    end do

    wInc = (wbar+be/R-be)*(1-w_size)/w_size
    ! call FindIndex(wInc,ybardist,mgrid,index,prob)
    ! AllMea = Mea*(1-w_size)
    ! AllMea(index-1) = AllMea(index-1)+w_size*prob
    ! AllMea(index) = AllMea(index)+w_size*(1-prob)

    ! construct inc_vec and mea_vec
    do indi = 1,mgrid
        do indj = 1,nz
            index = (indi-1)*nz+indj
            inc_vec(index) = ydist(indi,indj)
            mea_vec(index) = mea(indi)*pzvec(indj)
        end do
    end do
    mea_vec = mea_vec/sum(mea_vec)

    ! sort inc_vec and mea_vec by ascending order of inv_vec
    call piksr2(mgrid*nz,inc_vec,mea_vec)
    call FindIndex(wInc,inc_vec,mgrid*nz,index,prob)
    mea_vec = mea_vec*(1-w_size)
    mea_vec(index-1) = mea_vec(index-1)+w_size*prob
    mea_vec(index) = mea_vec(index)+w_size*(1-prob)

    ! calculate statistics
    perc = 0.99_dp
    stats(1) = FindPerc(perc,inc_vec,mea_vec,mgrid*nz)*100
    perc = 0.90_dp
    stats(2) = FindPerc(perc,inc_vec,mea_vec,mgrid*nz)*100
    stats(3) = FindGini(inc_vec,mea_vec,mgrid*nz)*100

end subroutine ss_distribution

subroutine t_step_distribution(Mea,b,bpr,p,ppr,MMea,stats)

    use parameter
    use global, only: zvec,pzvec,adist
    use routines

    implicit none
    real(dp), intent(in) :: b,p,bpr,ppr
    real(dp), dimension(mgrid), intent(in) :: Mea
    real(dp), dimension(mgrid), intent(out) :: MMea
    real(dp), dimension(3), intent(out) :: stats

    real(dp) :: w_size,be,phi,R,a,apr,wInc,perc,bepr
    real(dp), dimension(mgrid*nz) :: inc_vec,mea_vec
    real(dp), dimension(mgrid,nz) :: prodist,ydist
    integer, dimension(mgrid,nz) :: inddist
    integer :: indi, indj
    integer :: index, index_bar
    real(dp) :: prob, prob_bar

    w_size = wgt1

    be = w_size*b
    bepr = w_size*bpr
    phi = phifun(ppr,bepr)
    R = Rfun(be,bepr,ppr)

    do indi = 1,mgrid
        a = adist(indi)
        do indj = 1,nz
            apr = bbeta*((Afun(zVec(indj),ppr)+ppr)*phi/p+R*(1-phi))*a
            call FindIndex(apr,aDist,mgrid,index,prob)
            inddist(indi,indj) = index
            prodist(indi,indj) = prob
        end do
    end do

    ! new entrants
    apr = afun(zbar,ppr)+ppr+bepr
    call findindex(apr,adist,mgrid,index_bar,prob_bar)

    ! find next period distribution
    MMea=0.0

    do indi = 1,mgrid
        do indj = 1,nz
            index = indDist(indi,indj)
            prob = prodist(indi,indj)
            MMea(index-1) = MMea(index-1)+&
                surv_rate*Mea(indi)*prob*pzVec(indj)
            MMea(index) = MMea(index)+&
                surv_rate*Mea(indi)*(1-prob)*pzVec(indj)
        end do
    end do
    MMea(index_bar-1) = MMea(index_bar-1)+(1-surv_rate)*prob_bar
    MMea(index_bar) = MMea(index_bar)+(1-surv_rate)*(1-prob_bar)

    ! construct income distribution
    do indi = 1,mgrid
        a = adist(indi)
        ! ybardist(indi) = bbeta*(Afun(zbar,ppr)*phi/p+(R-1)*(1-phi))*a
        do indj = 1,nz
            ydist(indi,indj) = bbeta*(Afun(zVec(indj),ppr)*phi/p + &
                (R-1)*(1-phi))*a
        end do
    end do

    wInc = (wbar+bepr/R-bepr)*(1-w_size)/w_size
    ! call FindIndex(wInc,ybardist,mgrid,index,prob)
    ! AllMea = Mea*(1-w_size)
    ! AllMea(index-1) = AllMea(index-1)+w_size*prob
    ! AllMea(index) = AllMea(index)+w_size*(1-prob)

    ! construct inc_vec and mea_vec
    do indi = 1,mgrid
        do indj = 1,nz
            index = (indi-1)*nz+indj
            inc_vec(index) = ydist(indi,indj)
            mea_vec(index) = MMea(indi)*pzvec(indj)
        end do
    end do
    Mea_vec = Mea_vec/sum(Mea_vec)

    ! sort inc_vec and mea_vec by ascending order of inv_vec
    call piksr2(mgrid*nz,inc_vec,mea_vec)
    call FindIndex(wInc,inc_vec,mgrid*nz,index,prob)
    mea_vec = mea_vec*(1-w_size)
    mea_vec(index-1) = mea_vec(index-1)+w_size*prob
    mea_vec(index) = mea_vec(index)+w_size*(1-prob)

    ! calculate statistics
    perc = 0.99_dp
    stats(1) = FindPerc(perc,inc_vec,mea_vec,mgrid*nz)*100
    perc = 0.90_dp
    stats(2) = FindPerc(perc,inc_vec,mea_vec,mgrid*nz)*100
    stats(3) = FindGini(inc_vec,mea_vec,mgrid*nz)*100

end subroutine t_step_distribution