!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Subroutine SCMM-Hypo for Abaqus/Explicit
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Preprocessor definitions
!-----------------------------------------------------------------------
#include 'Definitions.f'
#ifndef SCMM_HYPO_EXPLICIT
#define SCMM_HYPO_EXPLICIT
!-----------------------------------------------------------------------
! Include files
!-----------------------------------------------------------------------
#include './Hypo.f'
!-----------------------------------------------------------------------
      subroutine vumat(
!-----------------------------------------------------------------------
! Read only (unmodifiable)variables -
!-----------------------------------------------------------------------
     +  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +  stepTime, totalTime, dt, cmname, coordMp, charLength,
     +  props, density, strainInc, relSpinInc,
     +  tempOld, stretchOld, defgradOld, fieldOld,
     +  stressOld, stateOld, enerInternOld, enerInelasOld,
     +  tempNew, stretchNew, defgradNew, fieldNew,
!-----------------------------------------------------------------------
! Write only (modifiable) variables -
!-----------------------------------------------------------------------
     +  stressNew, stateNew, enerInternNew, enerInelasNew )
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!      implicit none
      include 'vaba_param.inc'
!-----------------------------------------------------------------------
      integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
      real*8 stepTime, totalTime, dt
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +  charLength(nblock), strainInc(nblock,6),
     +  relSpinInc(nblock,3), tempOld(nblock),
     +  stretchOld(nblock,6),
     +  defgradOld(nblock,9),
     +  fieldOld(nblock,nfieldv), stressOld(nblock,6),
     +  stateOld(nblock,nstatev), enerInternOld(nblock),
     +  enerInelasOld(nblock), tempNew(nblock),
     +  stretchNew(nblock,6),
     +  defgradNew(nblock,9),
     +  fieldNew(nblock,nfieldv),
     +  stressNew(nblock,6), stateNew(nblock,nstatev),
     +  enerInternNew(nblock), enerInelasNew(nblock)
!-----------------------------------------------------------------------
      character*80 cmname
!-----------------------------------------------------------------------
!     Internal vumat variables
!-----------------------------------------------------------------------
      real*8 Fold(3,3)! Old Deformation gradient F=RU
      real*8 Fnew(3,3)! New Deformation gradient F=RU
      real*8 rr(3,3)! Rotation tensor R in polar decomposition F=RU
      real*8 rrt(3,3)! Transpose of rr
      real*8 uu(3,3)! Stretch tensor U in polar decomposition F=RU
      real*8 uui(3,3)! inverse of uu
      real*8 C11! Elastic coefficient
      real*8 C12! Elastic coefficient
      real*8 C44! Elastic coefficient
      real*8 xmat1(3,3), xmat2(3,3)! Tensors used for transformations
      real*8 sigsOld(nblock,6)! Old Stress tensor components, S11, S22, S33, S12, S23, S31 in global coordinate system
      real*8 sigsNew(nblock,6)! New Stress tensor components, S11, S22, S33, S12, S23, S31 in global coordinate system
      real*8 Dissipation(nblock)! The change in dissipated inelastic specific energy (sigma_ij*D^p_ij*dt=sum(tau(alpha)*dgamma(alpha)))
      real*8 zero, two, half
      integer km
      parameter (zero=0.d0, two=2.d0, half=5.d-1)
!-----------------------------------------------------------------------
!     This material subroutine is only for solid elements
!-----------------------------------------------------------------------
      if(ndir+nshr.ne.6)then
        call XPLB_ABQERR(-3,'This material subroutine is only for'//
     + ' solid elements',,,)
      endif
!-----------------------------------------------------------------------
!     Initial step
!-----------------------------------------------------------------------
      if(totaltime.eq.zero)then
!-----------------------------------------------------------------------
!       Read parameters from ABAQUS material card
!-----------------------------------------------------------------------
        C11        = props(1)! Elastic coefficient
        C12        = props(2)! Elastic coefficient
        C44        = props(3)! Elastic coefficient
!-----------------------------------------------------------------------
!       Assume small elastic deformation in the 
!       calculation of the initial elastic wave speeds
!-----------------------------------------------------------------------
        do km = 1, nblock
            stressNew(km,1) = C11*strainInc(km,1)+
     +                        C12*strainInc(km,2)+
     +                        C12*strainInc(km,3)
            stressNew(km,2) = C12*strainInc(km,1)+
     +                        C11*strainInc(km,2)+
     +                        C12*strainInc(km,3)
            stressNew(km,3) = C12*strainInc(km,1)+
     +                        C12*strainInc(km,2)+
     +                        C11*strainInc(km,3)
            stressNew(km,4) = two*C44*strainInc(km,4)
            stressNew(km,5) = two*C44*strainInc(km,5)
            stressNew(km,6) = two*C44*strainInc(km,6)
        enddo
        return
      endif
!-----------------------------------------------------------------------
!       Rotating the stress tensor to the Global coordinate frame
!       from rotated coordinate system used by Abaqus/Explicit
!-----------------------------------------------------------------------
      do km=1,nblock
!-----------------------------------------------------------------------
!       Old deformation gradient, F
!-----------------------------------------------------------------------
        Fold(1,1)   = defgradOld(km,1)
        Fold(2,2)   = defgradOld(km,2)
        Fold(3,3)   = defgradOld(km,3)
        Fold(1,2)   = defgradOld(km,4)
        Fold(2,3)   = defgradOld(km,5)
        Fold(3,1)   = defgradOld(km,6)
        Fold(2,1)   = defgradOld(km,7)
        Fold(3,2)   = defgradOld(km,8)
        Fold(1,3)   = defgradOld(km,9)
!-----------------------------------------------------------------------
!       Old stretch tensor, U
!-----------------------------------------------------------------------
        uu(1,1)     = stretchOld(km,1)
        uu(2,2)     = stretchOld(km,2)
        uu(3,3)     = stretchOld(km,3)
        uu(1,2)     = stretchOld(km,4)
        uu(2,3)     = stretchOld(km,5)
        uu(3,1)     = stretchOld(km,6)
        uu(2,1)     = stretchOld(km,4)
        uu(3,2)     = stretchOld(km,5)
        uu(1,3)     = stretchOld(km,6)
!-----------------------------------------------------------------------
!       Old Stress tensor in Rotated coordinate system used by Abaqus/Explicit
!-----------------------------------------------------------------------
        xmat1(1,1)  = stressOld(km,1)
        xmat1(2,2)  = stressOld(km,2)
        xmat1(3,3)  = stressOld(km,3)
        xmat1(1,2)  = stressOld(km,4)
        xmat1(2,3)  = stressOld(km,5)
        xmat1(3,1)  = stressOld(km,6)
        xmat1(2,1)  = stressOld(km,4)
        xmat1(3,2)  = stressOld(km,5)
        xmat1(1,3)  = stressOld(km,6)
!-----------------------------------------------------------------------
!       Transforming to Global coordinate system
!-----------------------------------------------------------------------
        call minv(uu,uui)
        call mmult(Fold,uui,rr)
        call mtransp(rr,rrt)
        call transform(xmat1,rr,rrt,xmat2)
        sigsOld(km,1) = xmat2(1,1)
        sigsOld(km,2) = xmat2(2,2)
        sigsOld(km,3) = xmat2(3,3)
        sigsOld(km,4) = xmat2(1,2)
        sigsOld(km,5) = xmat2(2,3)
        sigsOld(km,6) = xmat2(3,1)
      enddo
!-----------------------------------------------------------------------
!     Call the subroutine Hypo
!-----------------------------------------------------------------------
      call Hypo(sigsNew,stateNew,defgradNew,
     +          sigsOld,stateOld,defgradOld,dt,props,
     +          nblock,nstatev,nprops,Dissipation)
!-----------------------------------------------------------------------
!     Transforming the stress tensor from the global system to the Rotated coordinate system used in Abaqus/Explicit
!-----------------------------------------------------------------------
      do km=1,nblock
!-----------------------------------------------------------------------
!       New deformation gradient, F
!-----------------------------------------------------------------------
        Fnew(1,1)   = defgradNew(km,1)
        Fnew(2,2)   = defgradNew(km,2)
        Fnew(3,3)   = defgradNew(km,3)
        Fnew(1,2)   = defgradNew(km,4)
        Fnew(2,3)   = defgradNew(km,5)
        Fnew(3,1)   = defgradNew(km,6)
        Fnew(2,1)   = defgradNew(km,7)
        Fnew(3,2)   = defgradNew(km,8)
        Fnew(1,3)   = defgradNew(km,9)
!-----------------------------------------------------------------------
!       New stretch tensor, U
!-----------------------------------------------------------------------
        uu(1,1)     = stretchNew(km,1)
        uu(2,2)     = stretchNew(km,2)
        uu(3,3)     = stretchNew(km,3)
        uu(1,2)     = stretchNew(km,4)
        uu(2,3)     = stretchNew(km,5)
        uu(3,1)     = stretchNew(km,6)
        uu(2,1)     = stretchNew(km,4)
        uu(3,2)     = stretchNew(km,5)
        uu(1,3)     = stretchNew(km,6)
!-----------------------------------------------------------------------
!       New Stress tensor in Global coordinate system
!-----------------------------------------------------------------------
        xmat1(1,1)  = sigsNew(km,1)
        xmat1(2,2)  = sigsNew(km,2)
        xmat1(3,3)  = sigsNew(km,3)
        xmat1(1,2)  = sigsNew(km,4)
        xmat1(2,3)  = sigsNew(km,5)
        xmat1(3,1)  = sigsNew(km,6)
        xmat1(2,1)  = sigsNew(km,4)
        xmat1(3,2)  = sigsNew(km,5)
        xmat1(1,3)  = sigsNew(km,6)
!-----------------------------------------------------------------------
!       Transforming to rotated coordinate system used by Abaqus/Explicit
!-----------------------------------------------------------------------
        call minv(uu,uui)
        call mmult(Fnew,uui,rr)
        call mtransp(rr,rrt)
        call transform(xmat1,rrt,rr,xmat2)
        stressNew(km,1) = xmat2(1,1)
        stressNew(km,2) = xmat2(2,2)
        stressNew(km,3) = xmat2(3,3)
        stressNew(km,4) = xmat2(1,2)
        stressNew(km,5) = xmat2(2,3)
        stressNew(km,6) = xmat2(3,1)
      enddo
!-----------------------------------------------------------------------
!     Updating the specific internal energy
!-----------------------------------------------------------------------
      do km=1,nblock
        enerInternNew(km) = enerInternOld(km)+half*(
     +          (stressOld(km,1)+stressNew(km,1))*strainInc(km,1) +
     +          (stressOld(km,2)+stressNew(km,2))*strainInc(km,2) +
     +          (stressOld(km,3)+stressNew(km,3))*strainInc(km,3) +
     +      two*(stressOld(km,4)+stressNew(km,4))*strainInc(km,4) +
     +      two*(stressOld(km,5)+stressNew(km,5))*strainInc(km,5) +
     +      two*(stressOld(km,6)+stressNew(km,6))*strainInc(km,6))
     +      /density(km)
!-----------------------------------------------------------------------
!     Updating the dissipated inelastic specific energy
!-----------------------------------------------------------------------
        enerInelasNew(km) = enerInelasOld(km)+
     +                       Dissipation(km)/density(km)
      enddo
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      return
      end subroutine vumat
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------