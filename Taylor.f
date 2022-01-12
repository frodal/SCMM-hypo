!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Subroutine Taylor
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Subroutines should be inlined by the compiler
!-----------------------------------------------------------------------
!DIR$ ATTRIBUTES FORCEINLINE :: Taylor
!-----------------------------------------------------------------------
! Preprocessor definitions
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_TAYLOR
#define SCMM_HYPO_TAYLOR
#if SCMM_HYPO_DFLAG != 0
#define SCMM_HYPO_NSTATEV 30
#else
#define SCMM_HYPO_NSTATEV 28
#endif
!-----------------------------------------------------------------------
! Subroutine Taylor
!-----------------------------------------------------------------------
      subroutine Taylor(stressNew,stateNew,defgradNew,
     +               stressOld,stateOld,defgradOld,dt,props,
     +               nblock,nstatev,nprops,Dissipation)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nblock, nstatev, nprops
      real*8 dt
      real*8 props(nprops),
     +       stressOld(nblock,6), stateOld(nblock,nstatev),
     +       defgradNew(nblock,9), defgradOld(nblock,9),
     +       stressNew(nblock,6), stateNew(nblock,nstatev)
      real*8 Dissipation(nblock)! The change in dissipated inelastic 
!-----------------------------------------------------------------------
      integer i, j, k
      integer, parameter :: Nsdv = SCMM_HYPO_NSTATEV
      real*8 tempDissipation(nblock,8)
#if SCMM_HYPO_DFLAG != 0
      integer isActive ! Is the integration point active (0=deleted, 
!                                                         1=active)
#endif
!-----------------------------------------------------------------------
!     First step
!-----------------------------------------------------------------------
      if(stateold(1,13).lt.1.d-6)then ! First step
        do i=1,8
          do j=1,6
            do k=1,nblock
              stateOld(k,j+Nsdv+(i-1)*(Nsdv+6)) = stressOld(k,j)
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!     FC-Taylor homogenization for an 8 grain polycrystal
!-----------------------------------------------------------------------
      do i=1,8
!-----------------------------------------------------------------------
!       Call the subroutine Hypo
!-----------------------------------------------------------------------
      call Hypo(stateNew(:,1+Nsdv+(i-1)*(Nsdv+6):6+Nsdv+(i-1)*(Nsdv+6)),
     +          stateNew(:,1+(i-1)*(Nsdv+6):Nsdv+(i-1)*(Nsdv+6)),
     +          defgradNew,
     +          stateOld(:,1+Nsdv+(i-1)*(Nsdv+6):6+Nsdv+(i-1)*(Nsdv+6)),
     +          stateOld(:,1+(i-1)*(Nsdv+6):Nsdv+(i-1)*(Nsdv+6)),
     +          defgradOld,
     +          dt,props,nblock,Nsdv,nprops,tempDissipation(:,i))
      enddo
!-----------------------------------------------------------------------
!     FC-Taylor homogenization
!-----------------------------------------------------------------------
      stressNew = 0.d0
      do i=1,8
        do j=1,6
          do k=1,nblock
            stressNew(k,j) = stressNew(k,j) +
     +                       0.125*stateNew(k,j+Nsdv+(i-1)*(Nsdv+6))
          enddo
        enddo
      enddo
      Dissipation = 0.d0
      do i=1,8
        do k=1,nblock
          Dissipation(k) = Dissipation(k) + 0.125*tempDissipation(k,i)
        enddo
      enddo
!-----------------------------------------------------------------------
!     An element is deleted once the first FC-Taylor grain fails
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG != 0
      do k=1,nblock
        isActive = nint(min(stateNew(k,30),
     +                      stateNew(k,30+Nsdv+6),
     +                      stateNew(k,30+2*(Nsdv+6)),
     +                      stateNew(k,30+3*(Nsdv+6)),
     +                      stateNew(k,30+4*(Nsdv+6)),
     +                      stateNew(k,30+5*(Nsdv+6)),
     +                      stateNew(k,30+6*(Nsdv+6)),
     +                      stateNew(k,30+7*(Nsdv+6))))
        do i=1,8
          stateNew(k,30+(i-1)*(Nsdv+6)) = isActive
        enddo
      enddo
#endif
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      return
      end subroutine Taylor
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------