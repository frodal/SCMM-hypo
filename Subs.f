!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Subroutines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                         SUBROUTINE TRANSFORM
!-----------------------------------------------------------------------
! Rotates a matrix from one frame to another
! B = P.A.P^T
!-----------------------------------------------------------------------
      subroutine transform(a,xp,b)
!
      implicit none
!
      real*8, intent(in) :: a(3,3),xp(3,3)
      real*8, intent(out) :: b(3,3)
!     Local variables
      real*8 xpt(3,3),xmat(3,3)
!-----------------------------------------------------------------------
      call mtransp(xp,xpt)
      call mmult(a,xpt,xmat)
      call mmult(xp,xmat,b)
!
      return
      end subroutine transform
!
!-----------------------------------------------------------------------
!                         SUBROUTINE MINV
!-----------------------------------------------------------------------
! Computes the inverse of a (3x3) matrix AINV = A^-1
!-----------------------------------------------------------------------
      subroutine minv(a,ainv)
!
      implicit none
!
      real*8, intent(in) :: a(3,3)
      real*8, intent(out) :: ainv(3,3)
!     Local variables
      real*8 adj(3,3),det
      integer i,j
!-----------------------------------------------------------------------
! Compute the determinant of A
      call determ2(a,det)
! Compute the adjoint matrix of A
      adj(1,1) = a(2,2)*a(3,3)-a(3,2)*a(2,3)
      adj(1,2) = -(a(2,1)*a(3,3)-a(3,1)*a(2,3))
      adj(1,3) = a(2,1)*a(3,2)-a(3,1)*a(2,2)
      adj(2,1) = -(a(1,2)*a(3,3)-a(3,2)*a(1,3))
      adj(2,2) = a(1,1)*a(3,3)-a(3,1)*a(1,3)
      adj(2,3) = -(a(1,1)*a(3,2)-a(3,1)*a(1,2))
      adj(3,1) = a(1,2)*a(2,3)-a(2,2)*a(1,3)
      adj(3,2) = -(a(1,1)*a(2,3)-a(2,1)*a(1,3))
      adj(3,3) = a(1,1)*a(2,2)-a(2,1)*a(1,2)
! Compute the transpose of the adjoint matrix of A = inverse of A
      do j=1,3
        do i=1,3
          ainv(i,j) = adj(j,i)/det
        enddo
      enddo
!
      return
      end subroutine minv
!
!-----------------------------------------------------------------------
!                         SUBROUTINE MMULT
!-----------------------------------------------------------------------
! Computes the product of two (3x3) matrices AB = A.B
!-----------------------------------------------------------------------
      subroutine mmult(a,b,ab)
!
      implicit none
!
      real*8, intent(in) :: a(3,3),b(3,3)
      real*8, intent(out) :: ab(3,3)
!     Local variables
      integer i,j
!-----------------------------------------------------------------------
      do j=1,3
        do i=1,3
          ab(i,j) = a(i,1)*b(1,j)+a(i,2)*b(2,j)+a(i,3)*b(3,j)
        enddo
      enddo
!
      return
      end subroutine mmult
!
!-----------------------------------------------------------------------
!                         SUBROUTINE MTRANSP
!-----------------------------------------------------------------------
! Computes the transpose of a (3x3) matrix AT = A^T
!-----------------------------------------------------------------------
      subroutine mtransp(a,at)
!
      implicit none
!
      real*8, intent(in) :: a(3,3)
      real*8, intent(out) :: at(3,3)
!     Local variables
      integer i,j
!-----------------------------------------------------------------------
      do j=1,3
        do i=1,3
            at(i,j) = a(j,i)
        enddo
      enddo
!
      return
      end subroutine mtransp
!
!-----------------------------------------------------------------------
!                         SUBROUTINE MAT2VEC
!-----------------------------------------------------------------------
! Transforms a (3x3) symmetric matrix into a (6x1) vector
!-----------------------------------------------------------------------
      subroutine mat2vec(xmat,vec)
!
      implicit none
!
      real*8, intent(in) :: xmat(3,3)
      real*8, intent(out) :: vec(6)
!-----------------------------------------------------------------------
      vec(1) = xmat(1,1)
      vec(2) = xmat(2,2)
      vec(3) = xmat(3,3)
      vec(4) = xmat(1,2)
      vec(5) = xmat(2,3)
      vec(6) = xmat(3,1)
!
      return
      end subroutine mat2vec
!
!-----------------------------------------------------------------------
!                         SUBROUTINE VEC2MAT
!-----------------------------------------------------------------------
! Transforms a (6x1) vector into a (3x3) symmetric matrix
!-----------------------------------------------------------------------
      subroutine vec2mat(vec,xmat)
!
      implicit none
!
      real*8, intent(in) :: vec(6)
      real*8, intent(out) :: xmat(3,3)
!-----------------------------------------------------------------------
      xmat(1,1) = vec(1)
      xmat(2,2) = vec(2)
      xmat(3,3) = vec(3)
      xmat(1,2) = vec(4)
      xmat(2,1) = vec(4)
      xmat(2,3) = vec(5)
      xmat(3,2) = vec(5)
      xmat(1,3) = vec(6)
      xmat(3,1) = vec(6)
!
      return
      end subroutine vec2mat
!
!-----------------------------------------------------------------------
!                         SUBROUTINE DETERM
!-----------------------------------------------------------------------
! Computes the determinant of a (3x3) matrix   
!-----------------------------------------------------------------------
      subroutine determ2(a,det)
!
      implicit none
!
      real*8, intent(in) :: a(3,3)
      real*8, intent(out) :: det
!     Local variables
!-----------------------------------------------------------------------
      det = a(1,1)*a(2,2)*a(3,3)+
     +      a(2,1)*a(3,2)*a(1,3)+
     +      a(3,1)*a(1,2)*a(2,3)-
     +      a(3,1)*a(2,2)*a(1,3)-
     +      a(1,1)*a(3,2)*a(2,3)-
     +      a(2,1)*a(1,2)*a(3,3)
!
!      if (abs(det).lt.1.d-12) then
!        write(59,*)'Determinant is Zero!'
!        stop
!      endif
!
      return
      end subroutine determ2
!
!-----------------------------------------------------------------------
!                         SUBROUTINE SINC
!-----------------------------------------------------------------------
! Computes the strain and rotation increments from the deformation gradient F
! in the global coordinate system
!-----------------------------------------------------------------------
      subroutine sinc(Fold,Fnew,dt,epsinc,spininc)
!
      implicit none
!
      real*8, intent(in) :: Fnew(3,3), Fold(3,3), dt
      real*8, intent(out) :: epsinc(6), spininc(3)
!     Local variables
      real*8 L(3,3), Favginv(3,3), Favg(3,3), Fdot(3,3), Lt(3,3)
      real*8 Ddt(3,3), Wdt(3,3), half
      integer i, j
      parameter(half=5.d-1)
!-----------------------------------------------------------------------
      Fdot = (Fnew-Fold)/dt
      Favg = half*(Fnew+Fold)
      call minv(Favg,Favginv)
      call mmult(Fdot,Favginv,L)
      call mtransp(L,LT)
      Ddt       = half*dt*(L+LT)
      Wdt       = half*dt*(L-LT)
      call mat2vec(Ddt,epsinc)
      spininc(1)= Wdt(3,2)
      spininc(2)= Wdt(1,3)
      spininc(3)= Wdt(2,1)
!
      return
      end subroutine sinc
!
!-----------------------------------------------------------------------
!                         SUBROUTINE updateR
!-----------------------------------------------------------------------
! Computes the new rotation tensor based on lattice/elastic rotation increments
!-----------------------------------------------------------------------
      subroutine updateR(domega_e,R)
!
      implicit none
!
      real*8, intent(in) :: domega_e(3)
      real*8, intent(inout) :: R(3,3)
!     Local variables
      real*8 Ide(3,3), We(3,3), B1(3,3), B1i(3,3), B2(3,3), B(3,3)
      real*8 newR(3,3), zero, one, half
      parameter(zero=0.d0,one=1.d0,half=5.d-1)
!-----------------------------------------------------------------------
      ! Calculating B(i,l)=INV((I(i,k)-0.5*domega_e(i,k)))*(I(k,l)+0.5*domega_e(k,l))
      ! R(n+1)=B*R(n)
      We(1,1)  = zero
      We(1,2)  = -domega_e(3)
      We(1,3)  = domega_e(2)
      We(2,1)  = domega_e(3)
      We(2,2)  = zero
      We(2,3)  = -domega_e(1)
      We(3,1)  = -domega_e(2)
      We(3,2)  = domega_e(1)
      We(3,3)  = zero
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
!-----------------------------------------------------------------------
      B1       = Ide-half*We
      B2       = Ide+half*We
      call minv(B1,B1i)
      call mmult(B1i,B2,B)
      call mmult(B,R,newR)
      R        = newR
!-----------------------------------------------------------------------
      return
      end subroutine updateR
!
!-----------------------------------------------------------------------
!                         SUBROUTINE Voce
!-----------------------------------------------------------------------
! Update the critical resolved shear stresses/slip resistances
!-----------------------------------------------------------------------
      subroutine Voce(alpha,q,theta1,tau1,theta2,
     +                tau2,dgamma,gamma,tau_c)
!
      implicit none
!
      integer, intent(in) :: alpha
      real*8, intent(in) :: q(alpha,alpha),theta1,tau1,
     +                      theta2,tau2,dgamma(alpha),gamma
      real*8, intent(inout) :: tau_c(alpha)
!     Local variables
      integer a
!-----------------------------------------------------------------------
!     Voce
      do a=1,alpha
        tau_c(a)=tau_c(a)+(theta1*exp(-theta1*gamma/tau1)+
     +          theta2*exp(-theta2*gamma/tau2))*(q(a,1)*abs(dgamma(1))+
     +          q(a,2)*abs(dgamma(2))+q(a,3)*abs(dgamma(3))+
     +          q(a,4)*abs(dgamma(4))+q(a,5)*abs(dgamma(5))+
     +          q(a,6)*abs(dgamma(6))+q(a,7)*abs(dgamma(7))+
     +          q(a,8)*abs(dgamma(8))+q(a,9)*abs(dgamma(9))+
     +          q(a,10)*abs(dgamma(10))+q(a,11)*abs(dgamma(11))+
     +          q(a,12)*abs(dgamma(12)))
      enddo
!
      return
      end subroutine Voce
!
!-----------------------------------------------------------------------
!                         SUBROUTINE Kalidindi
!-----------------------------------------------------------------------
! Updates the critical resolved shear stresses/slip resistances
!-----------------------------------------------------------------------
      subroutine Kalidindi(alpha,q,h0,tau_s,am,dgamma,tau_c)
!
      implicit none
!
      integer, intent(in) :: alpha
      real*8, intent(in) :: q(alpha,alpha),h0,tau_s,am,
     +                      dgamma(alpha)
      real*8, intent(inout) :: tau_c(alpha)
!     Local variables
      real*8 dtau_c(alpha),one
      integer a
      parameter(one=1.d0)
!-----------------------------------------------------------------------
      do a=1,alpha
      ! Kalidindi et al.
      dtau_c(a)=q(a,1)*abs(dgamma(1))*h0*abs(one-tau_c(1)/tau_s)**am*
     .         sign(one,one-tau_c(1)/tau_s)+
     +         q(a,2)*abs(dgamma(2))*h0*abs(one-tau_c(2)/tau_s)**am*
     .         sign(one,one-tau_c(2)/tau_s)+
     +         q(a,3)*abs(dgamma(3))*h0*abs(one-tau_c(3)/tau_s)**am*
     .         sign(one,one-tau_c(3)/tau_s)+
     +         q(a,4)*abs(dgamma(4))*h0*abs(one-tau_c(4)/tau_s)**am*
     .         sign(one,one-tau_c(4)/tau_s)+
     +         q(a,5)*abs(dgamma(5))*h0*abs(one-tau_c(5)/tau_s)**am*
     .         sign(one,one-tau_c(5)/tau_s)+
     +         q(a,6)*abs(dgamma(6))*h0*abs(one-tau_c(6)/tau_s)**am*
     .         sign(one,one-tau_c(6)/tau_s)+
     +         q(a,7)*abs(dgamma(7))*h0*abs(one-tau_c(7)/tau_s)**am*
     .         sign(one,one-tau_c(7)/tau_s)+
     +         q(a,8)*abs(dgamma(8))*h0*abs(one-tau_c(8)/tau_s)**am*
     .         sign(one,one-tau_c(8)/tau_s)+
     +         q(a,9)*abs(dgamma(9))*h0*abs(one-tau_c(9)/tau_s)**am*
     .         sign(one,one-tau_c(9)/tau_s)+
     +         q(a,10)*abs(dgamma(10))*h0*abs(one-tau_c(10)/tau_s)**am*
     .         sign(one,one-tau_c(10)/tau_s)+
     +         q(a,11)*abs(dgamma(11))*h0*abs(one-tau_c(11)/tau_s)**am*
     .         sign(one,one-tau_c(11)/tau_s)+
     +         q(a,12)*abs(dgamma(12))*h0*abs(one-tau_c(12)/tau_s)**am*
     .         sign(one,one-tau_c(12)/tau_s)
      enddo
      tau_c = tau_c+dtau_c
!
      return
      end subroutine Kalidindi
!
!-----------------------------------------------------------------------
!                         SUBROUTINE EULER
!-----------------------------------------------------------------------
! Computes the Euler angles associated with the rotation matrix R
! in terms of the three Euler angles: phi1, PHI and phi2 (ang(1:3))
!-----------------------------------------------------------------------
      subroutine euler(R,ang)
!
      implicit none
!
      real*8, intent(in) :: R(3,3)
      real*8, intent(out) :: ang(3)
!     Local variables
      real*8 pi, zero, one, circ, small, halfCirc
      integer i
      parameter(pi=4.d0*atan(1.d0),zero=0.d0,one=1.d0,circ=360.d0,
     .          small=1.d-9,halfCirc=180.d0)
!-----------------------------------------------------------------------
      if (abs(abs(R(3,3))-one).gt.small) then
        ang(1) = atan2(R(1,3),-R(2,3))*halfCirc/pi
        ang(2) = acos(R(3,3))*halfCirc/pi
        ang(3) = atan2(R(3,1),R(3,2))*halfCirc/pi
      else
        ang(1) = atan2(R(2,1),R(1,1))*halfCirc/pi
        ang(2) = acos(R(3,3))*halfCirc/pi
        ang(3) = zero
      endif
!-----------------------------------------------------------------------
      if (ang(1).lt.zero) then
        ang(1) = ang(1)+circ
      endif
      if (ang(2).lt.zero) then
        ang(2) = ang(2)+circ
      endif
      if (ang(3).lt.zero) then
        ang(3) = ang(3)+circ
      endif
!
      return
      end subroutine euler