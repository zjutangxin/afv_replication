! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Main function for replicating "Financial Globalization, Inequality,
!       and the Rising Public Debt."
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

program debt_main

    use global
    use routines
    implicit none

    real(dp), dimension(nb,nb) :: DDebPol1,DDebPol2
    real(dp), dimension(nb,nb) :: DDebPol1fit,DDebPol2fit
    real(dp), dimension(7) :: coevec1,coevec2    
    real(dp), dimension(nb,nb) :: vv1wmx,vv2wmx,vv1emx,vv2emx
    real(dp), dimension(nb,nb) :: pprimx1,pprimx2
    real(dp), dimension(nt+1,nb,nb) :: Pri1_seqa, Pri2_seqa
    real(dp), dimension(nt+1,nb,nb) :: v1w_seqa,v2w_seqa
    real(dp), dimension(nt+1,nb,nb) :: v1e_seqa,v2e_seqa
    real(dp), dimension(nt,nb,nb) :: DebPol1_seqa, DebPol2_seqa
    real(dp), dimension(maxgrid,maxgrid) :: val1mx, val2mx
    real(dp) :: b1,b2,b1pr,b2pr,p1,p2,v1e,v1w,v2e,v2w
    real(dp), dimension(14) :: retval

    integer :: indb, indbp, indt, indb1, indb2, indz, indbp1, indbp2
    real(dp) :: fnorm, sumfnorm
    integer, dimension(1) :: maxind
    integer, dimension(maxgrid) :: maxindvec1,maxindvec2
    real(dp), dimension(2) :: debtGuess, fvec
    real(dp), dimension(nb,nb) :: GuessDebPol1, GuessDebPol2

    real(dp), dimension(mgrid) :: mea_ss
    real(dp) :: ssb, ssp
    real(dp), dimension(3) :: stats

    real(dp), parameter :: tol = 1e-15
    integer, parameter :: lwa = (2*(3*2+13))/2
    real(dp), dimension(lwa) :: wa
    integer :: info

    external NashSolution, find_ss
!=======================================================================!
!                          INITIALIZATION                               !
!=======================================================================!
    
    do indz = 1,nz
        zvec(indz) = zmin+(zmax-zmin)*(real(indz-1)/real(nz-1))        
    end do
    pzvec = 1.0_dp/real(nz)
    pzvec = pzvec/sum(pzvec)

    do indb = 1, nb
        bvec(indb) = bmin+(bmax-bmin)*(real(indb-1)/real(nb-1))
    end do

    do indbp = 1, maxGrid
        DebChoiceVec(indbp) = &
            bmin+(bmax-bmin)*(real(indbp-1)/real(maxGrid-1))
    end do    

    v1wMx = 0.0_dp
    v2wMx = 0.0_dp
    v1eMx = 0.0_dp
    v2eMx = 0.0_dp
    priMx1 = 0.0_dp
    priMx1 = 0.0_dp

    v1w_seqa(nt+1,:,:) = v1wMx
    v2w_seqa(nt+1,:,:) = v2wMx
    v1e_seqa(nt+1,:,:) = v1eMx
    v2e_seqa(nt+1,:,:) = v2eMx
    Pri1_seqa(nt+1,:,:) = PriMx1
    Pri2_seqa(nt+1,:,:) = PriMx2    

!=======================================================================!
!                          TERMINAL PERIOD                              !
!=======================================================================!

    indt = nt

    DebPol1 = 0.0_dp
    DebPol2 = 0.0_dp

    DebPol1_seqa(indt,:,:)=DebPol1
    DebPol2_seqa(indt,:,:)=DebPol2    

    do indb1 = 1, nb
        do indb2 = 1,nb
            b1 = bvec(indb1)
            b2 = bvec(indb2)
            b1pr = DebPol1(indb1,indb2)
            b2pr = DebPol2(indb1,indb2)
            call SolveSystem(indt,b1,b2,b1pr,b2pr,retval)
            p1  = RetVal(3)
            V1w = RetVal(5)
            V1e = RetVal(6)
            p2  = RetVal(9)
            V2w = RetVal(11)
            V2e = RetVal(12)
            vv1wMx(indb1,indb2)  = V1w
            vv2wMx(indb1,indb2)  = V2w
            vv1eMx(indb1,indb2)  = V1e
            vv2eMx(indb1,indb2)  = V2e
            PPriMx1(indb1,indb2) = p1
            PPriMx2(indb1,indb2) = p2

        end do
    end do

    v1wMx=vv1wMx
    v2wMx=vv2wMx
    v1eMx=vv1eMx
    v2eMx=vv2eMx
    PriMx1=PPriMx1
    PriMx2=PPriMx2

    v1w_seqa(indt,:,:)=vv1wMx
    v2w_seqa(indt,:,:)=vv2wMx
    v1e_seqa(indt,:,:)=vv1eMx
    v2e_seqa(indt,:,:)=vv2eMx
    Pri1_seqa(indt,:,:)=PPriMx1
    Pri2_seqa(indt,:,:)=PPriMx2

!=======================================================================!
!                       ITERATION STARTING AT T-1                       !
!=======================================================================!
    GuessDebPol1 = bbar
    GuessDebPol2 = bbar

    time: do indt = nt-1, 1, -1

        sumfnorm = 0.0_dp

        do indb1 = 1, nb
        do indb2 = 1, nb
            b1 = bvec(indb1)
            b2 = bvec(indb2)
            val1mx = 0.0_dp
            val2mx = 0.0_dp

            ! test parallel
            do indbp1 = 1,maxgrid
            do indbp2 = 1,maxgrid
                b1pr = DebChoiceVec(indbp1)
                b2pr = DebChoiceVec(indbp2)
                call SolveSystem(indt,b1,b2,b1pr,b2pr,RetVal)
                v1w = RetVal(5)
                v1e = RetVal(6)
                v2w = RetVal(11)
                v2e = RetVal(12)
                val1Mx(indbp1,indbp2) = v1w*wgt1+v1e*(1-wgt1)
                val2Mx(indbp1,indbp2) = v2w*wgt2+v2e*(1-wgt2)                
            enddo
            enddo
    
            ! best response functions
            do indbp2 = 1, MaxGrid
                MaxInd = MAXLOC(Val1Mx(:,indbp2))
                MaxIndVec1(indbp2) = MaxInd(1)
                RespFun1(indbp2)   = DebChoiceVec(MaxInd(1))
            end do
            
            do indbp1 = 1, MaxGrid
                MaxInd = MAXLOC(Val2Mx(indbp1,:))
                MaxIndVec2(indbp1) = MaxInd(1)
                RespFun2(indbp1)   = DebChoiceVec(MaxInd(1))
            end do

            ! solve the nash equilibrium
            debtGuess(1) = GuessDebPol1(indb1,indb2)
            debtGuess(2) = GuessDebPol2(indb1,indb2)
            fvec = 1000.0_dp
            call hybrd1(NashSolution,2,debtGuess,fvec,tol,info,wa,lwa)
            fnorm = sqrt(sum(fvec**2,1))
            SumFnorm = SumFnorm + Fnorm
            b1pr = debtGuess(1)
            b2pr = debtGuess(2)

            b1pr = min(b1pr,bmax)
            b2pr = min(b2pr,bmax)
            b1pr = max(b1pr,bmin)
            b2pr = max(b2pr,bmin)
            
            DDebPol1(indb1,indb2) = b1pr
            DDebPol2(indb1,indb2) = b2pr

        end do
        end do

        ! smooth the policy functions
        call PolynomialInt(bVec,Nb,bVec,Nb,DDebPol1,DDebPol1fit,CoeVec1)
        call PolynomialInt(bVec,Nb,bVec,Nb,DDebPol2,DDebPol2fit,CoeVec2)

        ! update guess
        DDebPol1 = DDebPol1fit
        DDebPol2 = DDebPol2fit
        GuessDebPol1 = DDebPol1
        GuessDebPol2 = DDebPol2

        write (*,'(a12,i6,3f17.11)') "POLICY ITER", nT-indt, &
        sum(abs(DDebPol1-DebPol1)),sum(abs(DDebPol2-DebPol2)),SumFnorm

        DebPol1 = DDebPol1
        DebPol2 = DDebPol2
        DebPol1_seqa(indt,:,:) = DebPol1
        DebPol2_seqa(indt,:,:) = DebPol2
        
        ! compute the equilibrium given Nash policies
        do indb1 = 1,nb
            do indb2 = 1,Nb
                b1 = bVec(indb1)
                b2 = bVec(indb2)
                b1pr = DebPol1(indb1,indb2)
                b2pr = DebPol2(indb1,indb2)
                call SolveSystem(indt,b1,b2,b1pr,b2pr,RetVal)
                p1  = RetVal(3)
                V1w = RetVal(5)
                V1e = RetVal(6)
                p2  = RetVal(9)
                V2w = RetVal(11)
                V2e = RetVal(12)

                vv1wMx(indb1,indb2) = V1w
                vv2wMx(indb1,indb2) = V2w
                vv1eMx(indb1,indb2) = V1e
                vv2eMx(indb1,indb2) = V2e
                PPriMx1(indb1,indb2) = p1
                PPriMx2(indb1,indb2) = p2
            end do
        end do
      
        v1wMx = vv1wMx
        v2wMx = vv2wMx
        v1eMx = vv1eMx
        v2eMx = vv2eMx
        PriMx1 = PPriMx1
        PriMx2 = PPriMx2
      
        v1w_seqa(indt,:,:) = v1wMx
        v2w_seqa(indt,:,:) = v2wMx
        v1e_seqa(indt,:,:) = v1eMx
        v2e_seqa(indt,:,:) = v2eMx
        Pri1_seqa(indt,:,:) = PriMx1
        Pri2_seqa(indt,:,:) = PriMx2        

        ! compute the implied steady state debt level when
        ! infinite horizon policies are approximated using
        ! t period problem
        debtGuess = bbar
        call hybrd1(find_ss,2,debtGuess,fvec,tol,info,wa,lwa)
        fnorm = sqrt(sum(fvec**2,1))
        write (*,101) 't = ', nt-indt, 'SS Debt 1 = ', debtGuess(1), &
            'SS Debt 2 = ', debtGuess(2), 'Error = ', fnorm
        101 format (A5, I3, A15, F10.6, A15, F10.6, A10, F12.6)            

        call SolveSystem(indt,debtGuess(1),debtGuess(2), &
            debtGuess(1),debtGuess(2),RetVal)
        write (*,'(A15,2F14.6)') 'Land Price = ', RetVal(3), RetVal(9)
        write (*,*) ''

    end do time

    open(1,file='./results/DebPol1.txt',form='formatted')
    do indb = 1,nb
        write(1,'(20ES14.6)') DebPol1(indb,:)
    end do
    close(1)

    open(1,file='./results/DebPol2.txt',form='formatted')
    do indb = 1,nb
        write(1,'(20ES14.6)') DebPol2(indb,:)
    end do
    close(1)    

!----------------------------------------------------------------!
!--------------COMPUTE TRANSITION INFINITE HORIZON---------------!
!----------------------------------------------------------------!
    ! ssb = debtGuess(1)
    ! ssp = retval(3)
    ! write (*,*) 'ssb = ', ssb, 'ssp = ', ssp

    ! test distribution
    ssb = 0.4616_dp
    ssp = 3.5062_dp

    call ss_distribution(ssb,ssp,Mea_ss,stats)

    write (*,*) 'Share top 1%,    Share top 10%,    Gini'
    write (*,'(3f14.6)') stats
    
end program debt_main