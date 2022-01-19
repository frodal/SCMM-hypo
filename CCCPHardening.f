!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Hardening subroutines for CCCP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Preprocessor definitions
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_CCCP_HARD
#define SCMM_HYPO_CCCP_HARD
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Subroutines
!-----------------------------------------------------------------------
! Subroutines should be inlined by the compiler
!-----------------------------------------------------------------------
!DIR$ ATTRIBUTES FORCEINLINE :: VoceMatrix, KalidindiMatrix, 
!DIR$& VoceCCCP, KalidindiCCCP
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                         SUBROUTINE VoceMatrix
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_KALIDINDI_ONLY
      subroutine VoceMatrix(q,theta1,tau1,theta2,tau2,gamma,hMatrix)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: q(alpha,alpha),theta1,tau1,theta2,tau2,gamma
      real*8, intent(out) :: hMatrix(alpha,alpha)
!     Local variables
      integer a,b
!-----
      do b=1,alpha
        do a=1,alpha
          hMatrix(a,b) = (theta1*exp(-theta1*gamma/tau1)+
     +          theta2*exp(-theta2*gamma/tau2))*q(a,b)
        enddo
      enddo
!
      return
      end subroutine VoceMatrix
#endif
!
!-----------------------------------------------------------------------
!                         SUBROUTINE KalidindiMatrix
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_VOCE_ONLY
      subroutine KalidindiMatrix(q,h0,tau_s,am,tau_c,hMatrix)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: q(alpha,alpha),h0,tau_s,am,tau_c(alpha)
      real*8, intent(out) :: hMatrix(alpha,alpha)
!     Local variables
      real*8, parameter :: one = 1.d0
      integer a,b
!-----
      do b=1,alpha
        do a=1,alpha
          hMatrix(a,b) = q(a,b)*h0*abs(one-tau_c(b)/tau_s)**am*
     .         sign(one,one-tau_c(b)/tau_s)
        enddo
      enddo
!
      return
      end subroutine KalidindiMatrix
#endif
!
!-----------------------------------------------------------------------
!                         SUBROUTINE VoceCCCP
!-----------------------------------------------------------------------
! Update the critical resolved shear stresses/slip resistances
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_KALIDINDI_ONLY
      subroutine VoceCCCP(q,theta1,tau1,theta2,
     +                tau2,dfdtau,dlambda,gamma,tau_c)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: q(alpha,alpha),theta1,tau1,
     +                      theta2,tau2,dfdtau(alpha),gamma,dlambda
      real*8, intent(inout) :: tau_c(alpha)
!     Local variables
      integer a
!-----
      do a=1,alpha
        ! Voce
        tau_c(a)=tau_c(a)+(theta1*exp(-theta1*gamma/tau1)+
     +          theta2*exp(-theta2*gamma/tau2))*(q(a,1)*abs(dfdtau(1))+
     +          q(a,2)*abs(dfdtau(2))+q(a,3)*abs(dfdtau(3))+
     +          q(a,4)*abs(dfdtau(4))+q(a,5)*abs(dfdtau(5))+
     +          q(a,6)*abs(dfdtau(6))+q(a,7)*abs(dfdtau(7))+
     +          q(a,8)*abs(dfdtau(8))+q(a,9)*abs(dfdtau(9))+
     +          q(a,10)*abs(dfdtau(10))+q(a,11)*abs(dfdtau(11))+
     +          q(a,12)*abs(dfdtau(12)))*dlambda
      enddo
!
      return
      end subroutine VoceCCCP
#endif
!
!-----------------------------------------------------------------------
!                         SUBROUTINE KalidindiCCCP
!-----------------------------------------------------------------------
! Updates the critical resolved shear stresses/slip resistances
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_VOCE_ONLY
      subroutine KalidindiCCCP(q,h0,tau_s,am,dfdtau,dlambda,tau_c)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: q(alpha,alpha),h0,tau_s,am,
     +                      dfdtau(alpha),dlambda
      real*8, intent(inout) :: tau_c(alpha)
!     Local variables
      real*8 dtau_c(alpha),one
      integer a
      parameter(one=1.d0)
!-----
      do a=1,alpha
        ! Kalidindi et al.
        dtau_c(a)=q(a,1)*abs(dfdtau(1))*h0*abs(one-tau_c(1)/tau_s)**am*
     .         sign(one,one-tau_c(1)/tau_s)+
     +         q(a,2)*abs(dfdtau(2))*h0*abs(one-tau_c(2)/tau_s)**am*
     .         sign(one,one-tau_c(2)/tau_s)+
     +         q(a,3)*abs(dfdtau(3))*h0*abs(one-tau_c(3)/tau_s)**am*
     .         sign(one,one-tau_c(3)/tau_s)+
     +         q(a,4)*abs(dfdtau(4))*h0*abs(one-tau_c(4)/tau_s)**am*
     .         sign(one,one-tau_c(4)/tau_s)+
     +         q(a,5)*abs(dfdtau(5))*h0*abs(one-tau_c(5)/tau_s)**am*
     .         sign(one,one-tau_c(5)/tau_s)+
     +         q(a,6)*abs(dfdtau(6))*h0*abs(one-tau_c(6)/tau_s)**am*
     .         sign(one,one-tau_c(6)/tau_s)+
     +         q(a,7)*abs(dfdtau(7))*h0*abs(one-tau_c(7)/tau_s)**am*
     .         sign(one,one-tau_c(7)/tau_s)+
     +         q(a,8)*abs(dfdtau(8))*h0*abs(one-tau_c(8)/tau_s)**am*
     .         sign(one,one-tau_c(8)/tau_s)+
     +         q(a,9)*abs(dfdtau(9))*h0*abs(one-tau_c(9)/tau_s)**am*
     .         sign(one,one-tau_c(9)/tau_s)+
     +         q(a,10)*abs(dfdtau(10))*h0*abs(one-tau_c(10)/tau_s)**am*
     .         sign(one,one-tau_c(10)/tau_s)+
     +         q(a,11)*abs(dfdtau(11))*h0*abs(one-tau_c(11)/tau_s)**am*
     .         sign(one,one-tau_c(11)/tau_s)+
     +         q(a,12)*abs(dfdtau(12))*h0*abs(one-tau_c(12)/tau_s)**am*
     .         sign(one,one-tau_c(12)/tau_s)
      enddo
      tau_c = tau_c+dtau_c*dlambda
!
      return
      end subroutine KalidindiCCCP
#endif
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------