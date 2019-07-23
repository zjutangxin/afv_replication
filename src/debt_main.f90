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

    use parameter
    implicit none

    real(dp), dimension(nbp) :: DebChoiceVec
    
    real(dp), dimension(nb,nb) :: DebPol1,DebPol2
    real(dp), dimension(nb,nb) :: vv1wmx,vv2wmx,vv1emx,vv2emx,pprimx1,pprimx2
    real(dp), dimension(nt+1,nb,nb) :: Pri1_seqa, Pri2_seqa
    real(dp), dimension(nt+1,nb,nb) :: v1w_seqa,v2w_seqa,v1e_seqa,v2e_seqa
    real(dp), dimension(nt,nb,nb) :: DebPol1_seqa, DebPol2_seqa    
    real(dp) :: b1, b2, b1pr, b2pr
    real(dp), dimension(14) :: retval

    integer :: indb, indbp, indt, indb1, indb2, indz

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
        DebChoiceVec(indb) = bmin+(bmax-bmin)*(real(indbp-1)/real(maxGrid-1))
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

    DebPol1_seqa(t,:,:)=DebPol1
    DebPol2_seqa(t,:,:)=DebPol2    

    do indb1 = 1, nb
        do indb2 = 1,nb
            b1 = bvec(indb1)
            b2 = bvec(indb2)
            b1pr = DebPol1(indb1,indb2)
            b2pr = DebPol2(indb1,indb2)
            call SolveSystem(t,b1,b2,b1pr,b2pr,retval)
            p1  = RetVal(3)
            V1w = RetVal(5)
            V1e = RetVal(6)
            p2  = RetVal(9)
            V2w = RetVal(11)
            V2e = RetVal(12)
            vv1wMx(i1,i2)  = V1w
            vv2wMx(i1,i2)  = V2w
            vv1eMx(i1,i2)  = V1e
            vv2eMx(i1,i2)  = V2e
            PPriMx1(i1,i2) = p1
            PPriMx2(i1,i2) = p2
        end do
    end do

    v1wMx=vv1wMx
    v2wMx=vv2wMx
    v1eMx=vv1eMx
    v2eMx=vv2eMx
    PriMx1=PPriMx1
	PriMx2=PPriMx2

    v1w_seqa(t,:,:)=vv1wMx
    v2w_seqa(t,:,:)=vv2wMx
    v1e_seqa(t,:,:)=vv1eMx
    v2e_seqa(t,:,:)=vv2eMx
    Pri1_seqa(t,:,:)=PPriMx1
    Pri2_seqa(t,:,:)=PPriMx2

end program debt_main