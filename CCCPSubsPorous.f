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
! Include files
!-----------------------------------------------------------------------
!DEC$ FREEFORM
#include './Dependencies/quartic.f90'
!DEC$ NOFREEFORM
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Subroutines
!-----------------------------------------------------------------------
! Subroutines should be inlined by the compiler
!-----------------------------------------------------------------------
!DIR$ ATTRIBUTES FORCEINLINE :: yieldfunction, yieldgradient, RMAP,
!DIR$& UpdateDamageHan, calcTauEff, findRoot
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                         SUBROUTINE yieldfunction
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine yieldfunction(tau_eff,tau_c,rho,m,f)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: tau_eff(alpha),tau_c(alpha),rho,m
      real*8, intent(out) :: f
!     Local variables
      real*8 temp
      real*8 zero,one
      parameter(zero=0.d0,one=1.d0)
      integer a
!-----
      temp = zero
      do a = 1,alpha
        temp = temp+exp((rho/m)*(abs(tau_eff(a))/(tau_c(a))-one))
      enddo
      f = (one/rho)*log(temp)
!
      return
      end subroutine yieldfunction
!
!-----------------------------------------------------------------------
!                         SUBROUTINE calcTauEff
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine calcTauEff(tau,sigma,VVF,aParam,q1,q2,tau_eff)
!
      use PolynomialRoots
      implicit none
!
      real*8, intent(in) :: tau(12),sigma(6),VVF,aParam,q1,q2
      real*8, intent(out) :: tau_eff(12)
!     Local variables
      integer alpha, a
      real*8 zero, one, half, three, oThree, precision
      parameter(alpha=12,zero=0.d0,one=1.d0,half=5.d-1,three=3.d0,
     .          oThree=one/three,precision=epsilon(zero))
      real*8 coeff(5), Sh, Seq2
      complex*16 z(4)
!-----------------------------------------------------------------------
      Seq2 = half*((sigma(1)-sigma(2))**2
     +      +(sigma(2)-sigma(3))**2
     +      +(sigma(3)-sigma(1))**2)
     +      +three*sigma(4)**2+three*sigma(5)**2
     +      +three*sigma(6)**2! Equivalent von Mises stress squared
      Sh = (sigma(1)+sigma(2)+sigma(3))*oThree ! hydrostatic stress
!-----------------------------------------------------------------------
      coeff(5) = 9.d0*q1*VVF*q2**8*Sh**8/358400000.d0
      coeff(4) = three*q1*VVF*q2**6*Sh**6/320000.d0
      coeff(3) = three*q1*VVF*q2**4*Sh**4/1600.d0
      coeff(1) = 2.d0*q1*VVF-q1**2*VVF**2-one
!-----------------------------------------------------------------------
      if(coeff(5).gt.precision)then
        do a=1, alpha
          coeff(2) = tau(a)**2+2.d0*aParam*VVF*Seq2/45.d0+
     +                       three*q1*VVF*q2**2*Sh**2/20.d0
          call QuarticRoots(coeff, z)
          call findRoot(tau, z, a, 4, tau_eff)
        enddo
      else
        do a=1, alpha
          tau_eff(a) = sqrt((tau(a)**2+2.d0*aParam*VVF*Seq2/45.d0+
     +    three*q1*VVF*q2**2*Sh**2/20.d0)/(one+
     +    q1**2*VVF**2-2.d0*q1*VVF))*sign(one,tau(a))
        enddo
      endif
!-----------------------------------------------------------------------
      return
      end subroutine calcTauEff
!
!-----------------------------------------------------------------------
!                         SUBROUTINE findRoot
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine findRoot(tau,z,a,Nroot,tau_eff)
!
      implicit none
!
      real*8, intent(in) :: tau(12)
      complex*16, intent(in) :: z(4)
      integer, intent(in) :: a, Nroot
      real*8, intent(out) :: tau_eff(12)
!     Local variables
      integer code, i
      real*8 zero, one, precision
      parameter(zero=0.d0,one=1.d0,precision=epsilon(zero))
!-----------------------------------------------------------------------
      code = 0
      do i=1,Nroot
        if((abs(aimag(z(i))).le.precision)
     .     .and.(dble(z(i)).gt.precision))then
          tau_eff(a) = one/sqrt(dble(z(i)))*sign(one,tau(a))
          code = code + 1
        endif
      enddo
      if(code.eq.0)then
        tau_eff(a) = tau(a)
#if defined SCMM_HYPO_STANDARD
        call STDB_ABQERR(-1,'A real value for tau_eff was not found'
     .                      ,,,)
#elif defined SCMM_HYPO_EXPLICIT
        call XPLB_ABQERR(-1,'A real value for tau_eff was not found'
     .                      ,,,)
#else
        write(*,*) 'Warning! A real value for tau_eff was not found'
#endif
      elseif(code.gt.1)then
#if defined SCMM_HYPO_STANDARD
        call STDB_ABQERR(-1,'More than one root found for tau_eff'
     .                      ,,,)
#elif defined SCMM_HYPO_EXPLICIT
        call XPLB_ABQERR(-1,'More than one root found for tau_eff'
     .                      ,,,)
#else
        write(*,*) 'Warning! More than one root found for tau_eff'
        write(*,*) "Number of roots", code
        write(*,*) " Roots: REAL PART   IMAGINARY PART"
        write(*,"(2ES20.12)") (dble(z(i)), aimag(z(i)), i=1,Nroot)
#endif
      endif
!-----------------------------------------------------------------------
      return
      end subroutine findRoot
!
!-----------------------------------------------------------------------
!                         SUBROUTINE yieldgradient
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine yieldgradient(tau,tau_eff,sigma,tau_c,VVF,rho,m,aParam,
     .                         q1,q2,S,dfdtau_eff,dfdtau_c,dfdsigma,
     .                         dfdsigmaskew,dfdVVF)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: tau(alpha),tau_eff(alpha),sigma(6),
     .                      tau_c(alpha),VVF,rho,m,aParam,q1,q2,
     .                      S(alpha,3,3)
      real*8, intent(out) :: dfdtau_eff(alpha),dfdtau_c(alpha),
     .                       dfdsigma(6),dfdsigmaskew(3),dfdVVF
!     Local variables
      real*8 dnom,temp,xmat(3,3),dgdtau_eff(12),dgdVVF(12),
     .       dgdsigma(12,3,3),sigXmat(3,3),Ide(3,3),Seq2,Sh
      real*8 zero,one,half,two,three,oThree,precision
      parameter(zero=0.d0,one=1.d0,half=5.d-1,two=2.d0,
     .          three=3.d0,oThree=one/three,precision=epsilon(zero))
      integer a,j,i
!-----------------------------------------------------------------------
      Seq2 = half*((sigma(1)-sigma(2))**2
     +      +(sigma(2)-sigma(3))**2
     +      +(sigma(3)-sigma(1))**2)
     +      +three*sigma(4)**2+three*sigma(5)**2
     +      +three*sigma(6)**2! Equivalent von Mises stress squared
      Sh = (sigma(1)+sigma(2)+sigma(3))*oThree ! hydrostatic stress
!-----------------------------------------------------------------------
      dnom   = zero
      do a = 1,alpha
        temp = exp((rho/m)*(abs(tau_eff(a))/(tau_c(a))-one))
        dnom = dnom+temp
        dfdtau_eff(a) = temp*sign(one,tau_eff(a))/(tau_c(a))
        dfdtau_c(a)   = -temp*abs(tau_eff(a))/((tau_c(a)**2))
      enddo
      dfdtau_eff = dfdtau_eff/(m*dnom)
      dfdtau_c   = dfdtau_c/(m*dnom)
      do a=1,alpha
        if(abs(tau_eff(a)).gt.precision)then
          dgdtau_eff(a) = -two*(tau(a)**2)/(tau_eff(a)**3)-
     .        4.d0*aParam*VVF*Seq2/(45.d0*tau_eff(a)**3)-
     . two*q1*VVF*(three*q2**2*Sh**2/(20.d0*tau_eff(a)**3)+
     .             three*q2**4*Sh**4/(800.d0*tau_eff(a)**5)+
     .             9.d0*q2**6*Sh**6/(320000.d0*tau_eff(a)**7)+
     .             9.d0*q2**8*Sh**8/(89600000.d0*tau_eff(a)**9))
          dgdVVF(a) = two*aParam*Seq2/(45.d0*tau_eff(a)**2)-
     .              two*q1**2*VVF+two*q1+
     . two*q1*(three*q2**2*Sh**2/(40.d0*tau_eff(a)**2)+
     .         three*q2**4*Sh**4/(3200.d0*tau_eff(a)**4)+
     .         three*q2**6*Sh**6/(640000.d0*tau_eff(a)**6)+
     .         9.d0*q2**8*Sh**8/(716800000.d0*tau_eff(a)**8))
        else
          dgdtau_eff(a) = one
          dgdVVF(a)     = zero
        endif
      enddo
      dfdVVF     = zero
      do a=1,alpha
        dfdVVF = dfdVVF - dfdtau_eff(a)*dgdVVF(a)/dgdtau_eff(a)
      enddo
!-----------------------------------------------------------------------
      Ide(1,1) = one
      Ide(1,2) = zero
      Ide(1,3) = zero
      Ide(2,1) = zero
      Ide(2,2) = one
      Ide(2,3) = zero
      Ide(3,1) = zero
      Ide(3,2) = zero
      Ide(3,3) = one
      call vec2mat(sigma,sigXmat)
!-----------------------------------------------------------------------
      do j=1,3
        do i=1,3
          do a=1,alpha
            if(abs(tau_eff(a)).gt.precision)then
              dgdsigma(a,i,j) = 
     . two*aParam*VVF*(sigXmat(i,j)-Sh*Ide(i,j))/(15.d0*tau_eff(a)**2)+
     . two*tau(a)*S(a,i,j)/tau_eff(a)**2+
     . (two/three)*q1*VVF*Ide(i,j)*(
     . three*q2**2*Sh/(20.d0*tau_eff(a)**2)+
     . three*q2**4*Sh**3/(800.d0*tau_eff(a)**4)+
     . 9.d0*q2**6*Sh**5/(320000.d0*tau_eff(a)**6)+
     . 9.d0*q2**8*Sh**7/(89600000.d0*tau_eff(a)**8))
            else
              dgdsigma(a,i,j) = zero
            endif
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      dfdsigmaskew = zero
      do a=1,alpha
        dfdsigmaskew(1) = dfdsigmaskew(1)-
     .                    half*dfdtau_eff(a)*(dgdsigma(a,3,2)-
     .                    dgdsigma(a,2,3))/dgdtau_eff(a)
        dfdsigmaskew(2) = dfdsigmaskew(2)-
     .                    half*dfdtau_eff(a)*(dgdsigma(a,1,3)-
     .                    dgdsigma(a,3,1))/dgdtau_eff(a)
        dfdsigmaskew(3) = dfdsigmaskew(3)-
     .                    half*dfdtau_eff(a)*(dgdsigma(a,2,1)-
     .                    dgdsigma(a,1,2))/dgdtau_eff(a)
      enddo
      xmat = zero
      do j=1,3
        do i=1,3
          do a=1,alpha
            xmat(i,j) = xmat(i,j)-
     .                  half*dfdtau_eff(a)*(dgdsigma(a,i,j)+
     .                  dgdsigma(a,j,i))/dgdtau_eff(a)
          enddo
        enddo
      enddo
      call mat2vec(xmat,dfdsigma)
!-----------------------------------------------------------------------
!
      return
      end subroutine yieldgradient
!
!-----------------------------------------------------------------------
!                         SUBROUTINE RMAP
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      subroutine RMAP(f,dfdtau_eff,dfdtau_c,dfdsigma,VVF,dfdVVF,
     .                sigma,C11,C12,C44,hMatrix,dlambda)
!
      implicit none
!
      integer, parameter :: alpha = 12
      real*8, intent(in) :: f,dfdtau_eff(alpha),dfdtau_c(alpha),
     .                      dfdsigma(6),VVF,dfdVVF,sigma(6),
     .                      C11,C12,C44,hMatrix(alpha,alpha)
      real*8, intent(out) :: dlambda
!     Local variables
      real*8 temp1,temp2,temp3,temp(6),Seq,Sh
      real*8 zero,half,one,three,four,oThree
      parameter(zero=0.d0,half=5.d-1,one=1.d0,three=3.d0,
     .          four=4.d0,oThree=1.d0/3.d0)
      integer a,b
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      temp1 = zero
      temp(1) = C11*(dfdsigma(1))*(one-VVF)
     +         +C12*(dfdsigma(2))*(one-VVF)
     +         +C12*(dfdsigma(3))*(one-VVF)
      temp(2) = C12*(dfdsigma(1))*(one-VVF)
     +         +C11*(dfdsigma(2))*(one-VVF)
     +         +C12*(dfdsigma(3))*(one-VVF)
      temp(3) = C12*(dfdsigma(1))*(one-VVF)
     +         +C12*(dfdsigma(2))*(one-VVF)
     +         +C11*(dfdsigma(3))*(one-VVF)
      temp(4) = four*C44*(dfdsigma(4))*(one-VVF)
      temp(5) = four*C44*(dfdsigma(5))*(one-VVF)
      temp(6) = four*C44*(dfdsigma(6))*(one-VVF)
      do a=1,6
        temp1 = temp1+temp(a)*dfdsigma(a)
      enddo
!-----------------------------------------------------------------------
      temp2 = zero
      do b=1,alpha
        do a=1,alpha
          temp2 = temp2+hMatrix(a,b)*abs(dfdtau_eff(b))*dfdtau_c(a)
        enddo
      enddo
!-----------------------------------------------------------------------
      temp3 = (one-VVF)**2*dfdVVF*(dfdsigma(1)+dfdsigma(2)+dfdsigma(3))
!-----------------------------------------------------------------------
      dlambda = f/(temp1-temp2-temp3)
!-----------------------------------------------------------------------
      return
      end subroutine RMAP
!
!-----------------------------------------------------------------------
!                         SUBROUTINE UpdateDamageHan
!-----------------------------------------------------------------------
! Updates the damage variable / void volume fraction
!-----------------------------------------------------------------------
      subroutine UpdateDamageHan(VVF,dfdsigma,dlambda)
!
      implicit none
!
      real*8, intent(inout) :: VVF
      real*8, intent(in) :: dfdsigma(6),dlambda
!     Local variables
      real*8 one,zero
      parameter(zero=0.d0,one=1.d0)
!-----------------------------------------------------------------------
      VVF = VVF + (one-VVF)**2*(dfdsigma(1)+dfdsigma(2)+
     .             dfdsigma(3))*dlambda
      VVF = max(VVF,zero)
!-----------------------------------------------------------------------
      return
      end subroutine UpdateDamageHan
!
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------