!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Subroutine SCMM-Hypo for Abaqus/Standard
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Subroutines should be inlined by the compiler
!-----------------------------------------------------------------------
!DIR$ ATTRIBUTES FORCEINLINE :: HypoImp3D, HypoImp2D
!-----------------------------------------------------------------------
! Preprocessor definitions
!-----------------------------------------------------------------------
#include 'Definitions.f'
#ifndef SCMM_HYPO_STANDARD
#define SCMM_HYPO_STANDARD
!-----------------------------------------------------------------------
! Include files
!-----------------------------------------------------------------------
#include './Taylor.f'
#include './Hypo.f'
#include './Subs.f'
!-----------------------------------------------------------------------
      subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     +                RPL,DDSDDT,DRPLDE,DRPLDT,
     +                STRAN,DSTRAN,TIME,dt,TEMP,DTEMP,PREDEF,DPRED,
     +                CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     +                DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     +                KSPT,JSTEP,KINC)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
      character*80 CMNAME
!-----------------------------------------------------------------------
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     +          DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     +          TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     +          DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!-----------------------------------------------------------------------
!     This material subroutine is only for ... elements
!-----------------------------------------------------------------------
#ifdef SCMM_HYPO_3D_ONLY
      if(ntens.ne.6)then
        call STDB_ABQERR(-3,'This material subroutine is only for'//
     + ' solid elements',,,)
      endif
      call HypoImp3D(STRESS, STATEV, DDSDDE, SSE, SPD,
     +                TIME, dt, NSTATV, PROPS, NPROPS,
     +                DROT, DFGRD0, DFGRD1)
#elif defined SCMM_HYPO_2D_ONLY
      if(ntens.ne.4)then
        call STDB_ABQERR(-3,'This material subroutine is only for'//
     + ' plane strain and axisymmetric elements',,,)
      endif
      call HypoImp2D(STRESS, STATEV, DDSDDE, SSE, SPD,
     +                TIME, dt, NSTATV, PROPS, NPROPS,
     +                DROT, DFGRD0, DFGRD1)
#else
      if(ntens.eq.6)then
        call HypoImp3D(STRESS, STATEV, DDSDDE, SSE, SPD,
     +                TIME, dt, NSTATV, PROPS, NPROPS,
     +                DROT, DFGRD0, DFGRD1)
      elseif(ntens.eq.4)then
        call HypoImp2D(STRESS, STATEV, DDSDDE, SSE, SPD,
     +                TIME, dt, NSTATV, PROPS, NPROPS,
     +                DROT, DFGRD0, DFGRD1)
      else
        call STDB_ABQERR(-3,'This material subroutine is only for'//
     + ' solid, plane strain and axisymmetric elements',,,)
      endif
#endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end subroutine UMAT
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                         SUBROUTINE HypoImp3D
!-----------------------------------------------------------------------
! A user material subroutine for Abaqus/Standard using solid elements
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_2D_ONLY
      subroutine HypoImp3D(STRESS, STATEV, DDSDDE, SSE, SPD,
     +                TIME, dt, NSTATV, PROPS, NPROPS,
     +                DROT, DFGRD0, DFGRD1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
      DIMENSION STRESS(6),STATEV(NSTATV),DDSDDE(6,6),
     +          TIME(2),PROPS(NPROPS),
     +          DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!-----------------------------------------------------------------------
      integer k,kk,i,j,flagCTO
      real*8 defgradOld(1,9),defgradNew(1,9)
      real*8 stateOld(1,NSTATV),stateNew(1,NSTATV)
      real*8 stressOld(1,6),stressNew(1,6)
      real*8 C11,C12,C44
      real*8 Dissipation(1)! The change in dissipated inelastic specific energy (sigma_ij*D^p_ij*dt=sum(tau(alpha)*dgamma(alpha)))
      real*8 DissipationTGT(12)
      real*8 SPRIME(6), TDROT(3,3)
      real*8 d(12,6),epsinc(6),spininc(3)
      real*8 L11(12),L12(12),L13(12)
      real*8 L21(12),L22(12),L23(12)
      real*8 L31(12),L32(12),L33(12)
      real*8 F(12,9),Fold(12,9)
      real*8 stateTGTold(12,NSTATV),stateTGTnew(12,NSTATV)
      real*8 stressTGTold(12,NSTATV),stressTGTnew(12,NSTATV)
      real*8 CTO(6,6),CTO2(6,6)
      real*8 pert,o2pert,zero,one,two,half
      parameter(pert=1e-6,o2pert=1.0/(2.0*pert),zero=0.d0,one=1.d0,
     .          two=2.d0,half=5.d-1)
!-----------------------------------------------------------------------
!     Packaging Deformation gradient
!-----------------------------------------------------------------------
      defgradOld(1,1) = DFGRD0(1,1)
      defgradOld(1,2) = DFGRD0(2,2)
      defgradOld(1,3) = DFGRD0(3,3)
      defgradOld(1,4) = DFGRD0(1,2)
      defgradOld(1,5) = DFGRD0(2,3)
      defgradOld(1,6) = DFGRD0(3,1)
      defgradOld(1,7) = DFGRD0(2,1)
      defgradOld(1,8) = DFGRD0(3,2)
      defgradOld(1,9) = DFGRD0(1,3)
!-----------------------------------------------------------------------
      defgradNew(1,1) = DFGRD1(1,1)
      defgradNew(1,2) = DFGRD1(2,2)
      defgradNew(1,3) = DFGRD1(3,3)
      defgradNew(1,4) = DFGRD1(1,2)
      defgradNew(1,5) = DFGRD1(2,3)
      defgradNew(1,6) = DFGRD1(3,1)
      defgradNew(1,7) = DFGRD1(2,1)
      defgradNew(1,8) = DFGRD1(3,2)
      defgradNew(1,9) = DFGRD1(1,3)
!-----------------------------------------------------------------------
!     Packaging history variables
!-----------------------------------------------------------------------
      do k=1,NSTATV
         stateOld(1,k) = STATEV(k)
      enddo
!-----------------------------------------------------------------------
!     Rotate the stress tensor to the correct coordinate system
!-----------------------------------------------------------------------
      call mtransp(DROT, TDROT)
      call ROTSIG(STRESS, TDROT, SPRIME, 1, 3, 3)
      STRESS = SPRIME
      stressOld(1,1) = STRESS(1)
      stressOld(1,2) = STRESS(2)
      stressOld(1,3) = STRESS(3)
      stressOld(1,4) = STRESS(4)
      stressOld(1,5) = STRESS(6)
      stressOld(1,6) = STRESS(5)
!-----------------------------------------------------------------------
!     Call the Subroutine
!-----------------------------------------------------------------------
#if SCMM_HYPO_MODEL == 3 || SCMM_HYPO_MODEL == 4
      call Taylor(stressNew, stateNew, defgradNew,
     +         stressOld, stateOld, defgradOld, dt, props,
     +         1, NSTATV, nprops, Dissipation)
#else
      call Hypo(stressNew, stateNew, defgradNew,
     +         stressOld, stateOld, defgradOld, dt, props,
     +         1, NSTATV, nprops, Dissipation)
#endif
!-----------------------------------------------------------------------
!     Update Consistent tangent operator
!-----------------------------------------------------------------------
      call sinc(DFGRD0, DFGRD1, dt, epsinc, spininc)
!-----------------------------------------------------------------------
      flagCTO = nint(props(17))
      if(flagCTO.eq.1)then
!-----------------------------------------------------------------------
!     Elastic tangent operator
!-----------------------------------------------------------------------
         C11 = props(1)
         C12 = props(2)
         C44 = props(3)
!-----------------------------------------------------------------------
         DDSDDE = zero
!-----------------------------------------------------------------------
         DDSDDE(1,1) = C11
         DDSDDE(1,2) = C12
         DDSDDE(1,3) = C12
!-----------------------------------------------------------------------
         DDSDDE(2,1) = C12
         DDSDDE(2,2) = C11
         DDSDDE(2,3) = C12
!-----------------------------------------------------------------------
         DDSDDE(3,1) = C12
         DDSDDE(3,2) = C12
         DDSDDE(3,3) = C11
!-----------------------------------------------------------------------
         DDSDDE(4,4) = C44
         DDSDDE(5,5) = C44
         DDSDDE(6,6) = C44
      else
!-----------------------------------------------------------------------
!     Consistent tangent operator
!-----------------------------------------------------------------------
         do i=1,6
            do k=1,12
               d(k,i) = epsinc(i)
            enddo
         enddo
!-----------------------------------------------------------------------
         kk = 1
         do i=1,6
            do k=1,2
               d(kk,i) = epsinc(i)+pert*(-one)**real(k)
               kk      = kk+1
            enddo
         enddo
!-----------------------------------------------------------------------
         do k=1,12
            L11(k) = d(k,1)
            L12(k) = d(k,4)-spininc(3)
            L13(k) = d(k,6)+spininc(2)
            L21(k) = d(k,4)+spininc(3)
            L22(k) = d(k,2)
            L23(k) = d(k,5)-spininc(1)
            L31(k) = d(k,6)-spininc(2)
            L32(k) = d(k,5)+spininc(1)
            L33(k) = d(k,3)
         enddo
!-----------------------------------------------------------------------
!     Update Consistent tangent operator
!-----------------------------------------------------------------------
         do k=1,12
            F(k,1) = DFGRD0(1,1)*(L11(k)+one)+DFGRD0(2,1)*L12(k)
     +           +DFGRD0(3,1)*L13(k)
            F(k,4) = DFGRD0(1,2)*(L11(k)+one)+DFGRD0(2,2)*L12(k)
     +           +DFGRD0(3,2)*L13(k)
            F(k,9) = DFGRD0(1,3)*(L11(k)+one)+DFGRD0(2,3)*L12(k)
     +           +DFGRD0(3,3)*L13(k)
            F(k,7) = DFGRD0(2,1)*(L22(k)+one)+DFGRD0(1,1)*L21(k)
     +           +DFGRD0(3,1)*L23(k)
            F(k,2) = DFGRD0(2,2)*(L22(k)+one)+DFGRD0(1,2)*L21(k)
     +           +DFGRD0(3,2)*L23(k)
            F(k,5) = DFGRD0(2,3)*(L22(k)+one)+DFGRD0(1,3)*L21(k)
     +           +DFGRD0(3,3)*L23(k)
            F(k,6) = DFGRD0(3,1)*(L33(k)+one)+DFGRD0(1,1)*L31(k)
     +           +DFGRD0(2,1)*L32(k)
            F(k,8) = DFGRD0(3,2)*(L33(k)+one)+DFGRD0(1,2)*L31(k)
     +           +DFGRD0(2,2)*L32(k)
            F(k,3) = DFGRD0(3,3)*(L33(k)+one)+DFGRD0(1,3)*L31(k)
     +           +DFGRD0(2,3)*L32(k)
         enddo
!-----------------------------------------------------------------------
         do k=1,12
            Fold(k,1) = DFGRD0(1,1)
            Fold(k,2) = DFGRD0(2,2)
            Fold(k,3) = DFGRD0(3,3)
            Fold(k,4) = DFGRD0(1,2)
            Fold(k,5) = DFGRD0(2,3)
            Fold(k,6) = DFGRD0(3,1)
            Fold(k,7) = DFGRD0(2,1)
            Fold(k,8) = DFGRD0(3,2)
            Fold(k,9) = DFGRD0(1,3)
         enddo
!-----------------------------------------------------------------------
         do i=1,NSTATV
            do k=1,12
               stateTGTold(k,i) = STATEV(i)
            enddo
         enddo
!-----------------------------------------------------------------------
          do k=1,12
              stressTGTold(k,1) = STRESS(1)
              stressTGTold(k,2) = STRESS(2)
              stressTGTold(k,3) = STRESS(3)
              stressTGTold(k,4) = STRESS(4)
              stressTGTold(k,5) = STRESS(6)
              stressTGTold(k,6) = STRESS(5)
          enddo
!-----------------------------------------------------------------------
!        Calculating stress state based on perturbation
!-----------------------------------------------------------------------
#if SCMM_HYPO_MODEL == 3 || SCMM_HYPO_MODEL == 4
         call Taylor(stressTGTnew, stateTGTnew, F,
     +         stressTGTold, stateTGTold, Fold, dt, props,
     +         12, NSTATV, nprops, DissipationTGT)
#else
         call Hypo(stressTGTnew, stateTGTnew, F,
     +         stressTGTold, stateTGTold, Fold, dt, props,
     +         12, NSTATV, nprops, DissipationTGT)
#endif
!-----------------------------------------------------------------------
         kk = 0
         do i=1,12,2
            kk = kk+1
            do j=1,6
               CTO(kk,j) = (stressTGTnew(i+1,j)-stressTGTnew(i,j))
     .                     *o2pert
            enddo
         enddo
!-----------------------------------------------------------------------
         do i=1,6
            do j=1,6
               CTO2(j,i) = CTO(i,j)
            enddo
         enddo
!-----------------------------------------------------------------------
         do j=1,4
            do i=1,4
               DDSDDE(i,j) = CTO2(i,j)
            enddo
         enddo
         do j=1,4
            DDSDDE(5,j) = CTO2(6,j)
            DDSDDE(6,j) = CTO2(5,j)
            DDSDDE(j,5) = CTO2(j,6)
            DDSDDE(j,6) = CTO2(j,5)
         enddo
         DDSDDE(5,5) = CTO2(6,6)
         DDSDDE(5,6) = CTO2(6,5)
         DDSDDE(6,5) = CTO2(5,6)
         DDSDDE(6,6) = CTO2(5,5)
      endif
!-----------------------------------------------------------------------
!     Update history variables
!-----------------------------------------------------------------------
      do k=1,NSTATV
         STATEV(k) = stateNew(1,k)
      enddo
!-----------------------------------------------------------------------
!     Update stress
!-----------------------------------------------------------------------
      STRESS(1) = stressNew(1,1)
      STRESS(2) = stressNew(1,2)
      STRESS(3) = stressNew(1,3)
      STRESS(4) = stressNew(1,4)
      STRESS(6) = stressNew(1,5)
      STRESS(5) = stressNew(1,6)
!-----------------------------------------------------------------------
!     Update Energies
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Updating the dissipated inelastic specific energy
!-----------------------------------------------------------------------
      SPD = SPD+Dissipation(1)
!-----------------------------------------------------------------------
!     Updating the specific elastic internal energy
!-----------------------------------------------------------------------
      SSE = SSE+(half*(
     +          (stressOld(1,1)+stressNew(1,1))*epsinc(1) +
     +          (stressOld(1,2)+stressNew(1,2))*epsinc(2) +
     +          (stressOld(1,3)+stressNew(1,3))*epsinc(3) +
     +      two*(stressOld(1,4)+stressNew(1,4))*epsinc(4) +
     +      two*(stressOld(1,5)+stressNew(1,5))*epsinc(5) +
     +      two*(stressOld(1,6)+stressNew(1,6))*epsinc(6))
     +       -Dissipation(1))
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end subroutine HypoImp3D
#endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                         SUBROUTINE HypoImp2D
!-----------------------------------------------------------------------
! A user material subroutine for Abaqus/Standard using plane strain and 
! axisymmetric elements
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_3D_ONLY
      subroutine HypoImp2D(STRESS, STATEV, DDSDDE, SSE, SPD,
     +                TIME, dt, NSTATV, PROPS, NPROPS,
     +                DROT, DFGRD0, DFGRD1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
      DIMENSION STRESS(4),STATEV(NSTATV),DDSDDE(4,4),
     +          TIME(2),PROPS(NPROPS),
     +          DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!-----------------------------------------------------------------------
      integer k,kk,i,j,flagCTO
      real*8 defgradOld(1,9),defgradNew(1,9)
      real*8 stateOld(1,NSTATV),stateNew(1,NSTATV)
      real*8 stressOld(1,6),stressNew(1,6)
      real*8 C11,C12,C44
      real*8 Dissipation(1)! The change in dissipated inelastic specific energy (sigma_ij*D^p_ij*dt=sum(tau(alpha)*dgamma(alpha)))
      real*8 DissipationTGT(8)
      real*8 SPRIME(4), TDROT(3,3)
      real*8 d(8,4),epsinc(6),spininc(3)
      real*8 L11(8),L12(8)
      real*8 L21(8),L22(8)
      real*8 L33(8)
      real*8 F(8,9),Fold(8,9)
      real*8 stateTGTold(8,NSTATV),stateTGTnew(8,NSTATV)
      real*8 stressTGTold(8,NSTATV),stressTGTnew(8,NSTATV)
      real*8 pert,o2pert,zero,one,two,half
      parameter(pert=1e-6,o2pert=1.0/(2.0*pert),zero=0.d0,one=1.d0,
     .          two=2.d0,half=5.d-1)
!-----------------------------------------------------------------------
!     Initial step
!-----------------------------------------------------------------------
      if(time(1).eq.zero.and.time(2).eq.zero)then
        call STDB_ABQERR(-1,'This material subroutine is only for'//
     + ' plane strain and axisymmetric elements with'//
     + ' certain crystallographic orientations.'//
     + ' Do not use it with plane stress elements!',,,)
      endif
!-----------------------------------------------------------------------
!     Packaging Deformation gradient
!-----------------------------------------------------------------------
      defgradOld(1,1) = DFGRD0(1,1)
      defgradOld(1,2) = DFGRD0(2,2)
      defgradOld(1,3) = DFGRD0(3,3)
      defgradOld(1,4) = DFGRD0(1,2)
      defgradOld(1,5) = DFGRD0(2,3)
      defgradOld(1,6) = DFGRD0(3,1)
      defgradOld(1,7) = DFGRD0(2,1)
      defgradOld(1,8) = DFGRD0(3,2)
      defgradOld(1,9) = DFGRD0(1,3)
!-----------------------------------------------------------------------
      defgradNew(1,1) = DFGRD1(1,1)
      defgradNew(1,2) = DFGRD1(2,2)
      defgradNew(1,3) = DFGRD1(3,3)
      defgradNew(1,4) = DFGRD1(1,2)
      defgradNew(1,5) = DFGRD1(2,3)
      defgradNew(1,6) = DFGRD1(3,1)
      defgradNew(1,7) = DFGRD1(2,1)
      defgradNew(1,8) = DFGRD1(3,2)
      defgradNew(1,9) = DFGRD1(1,3)
!-----------------------------------------------------------------------
!     Packaging history variables
!-----------------------------------------------------------------------
      do k=1,NSTATV
         stateOld(1,k) = STATEV(k)
      enddo
!-----------------------------------------------------------------------
!     Rotate the stress tensor to the correct coordinate system
!-----------------------------------------------------------------------
      call mtransp(DROT, TDROT)
      call ROTSIG(STRESS, TDROT, SPRIME, 1, 3, 1)
      STRESS = SPRIME
      stressOld(1,1) = STRESS(1)
      stressOld(1,2) = STRESS(2)
      stressOld(1,3) = STRESS(3)
      stressOld(1,4) = STRESS(4)
      stressOld(1,5) = zero
      stressOld(1,6) = zero
!-----------------------------------------------------------------------
!     Call the Subroutine
!-----------------------------------------------------------------------
#if SCMM_HYPO_MODEL == 3 || SCMM_HYPO_MODEL == 4
      call Taylor(stressNew, stateNew, defgradNew,
     +         stressOld, stateOld, defgradOld, dt, props,
     +         1, NSTATV, nprops, Dissipation)
#else
      call Hypo(stressNew, stateNew, defgradNew,
     +         stressOld, stateOld, defgradOld, dt, props,
     +         1, NSTATV, nprops, Dissipation)
#endif
!-----------------------------------------------------------------------
!     Update Consistent tangent operator
!-----------------------------------------------------------------------
      call sinc(DFGRD0, DFGRD1, dt, epsinc, spininc)
!-----------------------------------------------------------------------
      flagCTO = nint(props(17))
      if(flagCTO.eq.1)then
!-----------------------------------------------------------------------
!     Elastic tangent operator
!-----------------------------------------------------------------------
         C11 = props(1)
         C12 = props(2)
         C44 = props(3)
!-----------------------------------------------------------------------
         DDSDDE = zero
!-----------------------------------------------------------------------
         DDSDDE(1,1) = C11
         DDSDDE(1,2) = C12
         DDSDDE(1,3) = C12
!-----------------------------------------------------------------------
         DDSDDE(2,1) = C12
         DDSDDE(2,2) = C11
         DDSDDE(2,3) = C12
!-----------------------------------------------------------------------
         DDSDDE(3,1) = C12
         DDSDDE(3,2) = C12
         DDSDDE(3,3) = C11
!-----------------------------------------------------------------------
         DDSDDE(4,4) = C44
      else
!-----------------------------------------------------------------------
!     Consistent tangent operator
!-----------------------------------------------------------------------
         do i=1,4
            do k=1,8
               d(k,i) = epsinc(i)
            enddo
         enddo
!-----------------------------------------------------------------------
         kk = 1
         do i=1,4
            do k=1,2
               d(kk,i) = epsinc(i)+pert*(-one)**real(k)
               kk      = kk+1
            enddo
         enddo
!-----------------------------------------------------------------------
         do k=1,8
            L11(k) = d(k,1)
            L12(k) = d(k,4)-spininc(3)
            L21(k) = d(k,4)+spininc(3)
            L22(k) = d(k,2)
            L33(k) = d(k,3)
         enddo
!-----------------------------------------------------------------------
!     Update Consistent tangent operator
!-----------------------------------------------------------------------
         do k=1,8
            F(k,1) = DFGRD0(1,1)*(L11(k)+one)+DFGRD0(2,1)*L12(k)
            F(k,4) = DFGRD0(1,2)*(L11(k)+one)+DFGRD0(2,2)*L12(k)
            F(k,9) = zero
            F(k,7) = DFGRD0(2,1)*(L22(k)+one)+DFGRD0(1,1)*L21(k)
            F(k,2) = DFGRD0(2,2)*(L22(k)+one)+DFGRD0(1,2)*L21(k)
            F(k,5) = zero
            F(k,6) = zero
            F(k,8) = zero
            F(k,3) = DFGRD0(3,3)*(L33(k)+one)
         enddo
!-----------------------------------------------------------------------
         do k=1,8
            Fold(k,1) = DFGRD0(1,1)
            Fold(k,2) = DFGRD0(2,2)
            Fold(k,3) = DFGRD0(3,3)
            Fold(k,4) = DFGRD0(1,2)
            Fold(k,5) = zero
            Fold(k,6) = zero
            Fold(k,7) = DFGRD0(2,1)
            Fold(k,8) = zero
            Fold(k,9) = zero
         enddo
!-----------------------------------------------------------------------
         do i=1,NSTATV
            do k=1,8
               stateTGTold(k,i) = STATEV(i)
            enddo
         enddo
!-----------------------------------------------------------------------
          do k=1,8
              stressTGTold(k,1) = STRESS(1)
              stressTGTold(k,2) = STRESS(2)
              stressTGTold(k,3) = STRESS(3)
              stressTGTold(k,4) = STRESS(4)
              stressTGTold(k,5) = zero
              stressTGTold(k,6) = zero
          enddo
!-----------------------------------------------------------------------
!        Calculating stress state based on perturbation
!-----------------------------------------------------------------------
#if SCMM_HYPO_MODEL == 3 || SCMM_HYPO_MODEL == 4
         call Taylor(stressTGTnew, stateTGTnew, F,
     +         stressTGTold, stateTGTold, Fold, dt, props,
     +         8, NSTATV, nprops, DissipationTGT)
#else
         call Hypo(stressTGTnew, stateTGTnew, F,
     +         stressTGTold, stateTGTold, Fold, dt, props,
     +         8, NSTATV, nprops, DissipationTGT)
#endif
!-----------------------------------------------------------------------
         kk = 0
         do i=1,8,2
            kk = kk+1
            do j=1,4
                DDSDDE(j,kk) = (stressTGTnew(i+1,j)-stressTGTnew(i,j))
     .                     *o2pert
            enddo
         enddo
      endif
!-----------------------------------------------------------------------
!     Update history variables
!-----------------------------------------------------------------------
      do k=1,NSTATV
         STATEV(k) = stateNew(1,k)
      enddo
!-----------------------------------------------------------------------
!     Update stress
!-----------------------------------------------------------------------
      STRESS(1) = stressNew(1,1)
      STRESS(2) = stressNew(1,2)
      STRESS(3) = stressNew(1,3)
      STRESS(4) = stressNew(1,4)
!-----------------------------------------------------------------------
!     Update Energies
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Updating the dissipated inelastic specific energy
!-----------------------------------------------------------------------
      SPD = SPD+Dissipation(1)
!-----------------------------------------------------------------------
!     Updating the specific elastic internal energy
!-----------------------------------------------------------------------
      SSE = SSE+(half*(
     +          (stressOld(1,1)+stressNew(1,1))*epsinc(1) +
     +          (stressOld(1,2)+stressNew(1,2))*epsinc(2) +
     +          (stressOld(1,3)+stressNew(1,3))*epsinc(3) +
     +      two*(stressOld(1,4)+stressNew(1,4))*epsinc(4))
     +       -Dissipation(1))
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end subroutine HypoImp2D
#endif
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------