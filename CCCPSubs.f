!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Helper subroutines for CCCP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Preprocessor definitions
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_CCCP_SUBS
#define SCMM_HYPO_CCCP_SUBS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Subroutines
!-----------------------------------------------------------------------
! Subroutines should be inlined by the compiler
!-----------------------------------------------------------------------
!DIR$ ATTRIBUTES FORCEINLINE :: yieldfunction, yieldgradient, RMAP
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                         SUBROUTINE yieldfunction
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine yieldfunction(tau,tau_c,rho,m,f)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: tau(alpha),tau_c(alpha),rho,m
      real*8, intent(out) :: f
!     Local variables
      real*8 temp
      real*8 zero,one
      parameter(zero=0.d0,one=1.d0)
      integer a
!-----
      temp = zero
      do a = 1,alpha
        temp = temp+exp((rho/m)*(abs(tau(a))/tau_c(a)-one))
      enddo
      f = (one/rho)*log(temp)
!
      return
      end subroutine yieldfunction
!
!-----------------------------------------------------------------------
!                         SUBROUTINE yieldgradient
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine yieldgradient(tau,tau_c,rho,m,
     .                         S,dfdtau,dfdtau_c,dfdsigma)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: tau(alpha),tau_c(alpha),rho,m,S(alpha,3,3)
      real*8, intent(out) :: dfdtau(alpha),dfdtau_c(alpha),dfdsigma(6)
!     Local variables
      real*8 temp,xmat(3,3)
      real*8 zero,one,half
      parameter(zero=0.d0,one=1.d0,half=5.d-1)
      integer a,j,i
!-----
      temp = zero
      do a = 1,alpha
        temp = temp+exp((rho/m)*(abs(tau(a))/tau_c(a)-one))
        dfdtau(a) = (sign(one,tau(a))/tau_c(a))*
     .               exp((rho/m)*(abs(tau(a))/tau_c(a)-one))
        dfdtau_c(a) = -(abs(tau(a))/(tau_c(a))**2)*
     .                  exp((rho/m)*(abs(tau(a))/tau_c(a)-one))
      enddo
      dfdtau    = dfdtau/(m*temp)
      dfdtau_c  = dfdtau_c/(m*temp)
      xmat      = zero
      do j=1,3
        do i=1,3
          do a=1,alpha
            xmat(i,j)=xmat(i,j)+half*dfdtau(a)*(S(a,i,j)+S(a,j,i))
          enddo
        enddo
      enddo
      call mat2vec(xmat,dfdsigma)
!
      return
      end subroutine yieldgradient
!
!-----------------------------------------------------------------------
!                         SUBROUTINE RMAP
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine RMAP(f,dfdtau,dfdtau_c,dfdsigma,C11,C12,C44,hMatrix,
     .                dlambda)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: f,dfdtau(alpha),dfdtau_c(alpha),
     .                      dfdsigma(6),C11,C12,C44,
     .                      hMatrix(alpha,alpha)
      real*8, intent(out) :: dlambda
!     Local variables
      real*8 temp1,temp2,temp(6)
      real*8 zero,four
      parameter(zero=0.d0,four=4.d0)
      integer a,b
!-----
      temp1 = zero
      temp2 = zero
      temp(1) = C11*(dfdsigma(1))
     +         +C12*(dfdsigma(2))
     +         +C12*(dfdsigma(3))
      temp(2) = C12*(dfdsigma(1))
     +         +C11*(dfdsigma(2))
     +         +C12*(dfdsigma(3))
      temp(3) = C12*(dfdsigma(1))
     +         +C12*(dfdsigma(2))
     +         +C11*(dfdsigma(3))
      temp(4) = four*C44*(dfdsigma(4))
      temp(5) = four*C44*(dfdsigma(5))
      temp(6) = four*C44*(dfdsigma(6))
      do a=1,6
        temp1 = temp1+temp(a)*dfdsigma(a)
      enddo
      do b=1,alpha
        do a=1,alpha
          temp2 = temp2+hMatrix(a,b)*abs(dfdtau(b))*dfdtau_c(a)
        enddo
      enddo
      dlambda = f/(temp1-temp2)
!
      return
      end subroutine RMAP
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------