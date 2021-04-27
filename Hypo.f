!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Subroutine Hypo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Subroutines should be inlined by the compiler
!-----------------------------------------------------------------------
!DIR$ ATTRIBUTES FORCEINLINE :: Hypo
!-----------------------------------------------------------------------
! Preprocessor definitions
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_HYPO
#define SCMM_HYPO_HYPO
!-----------------------------------------------------------------------
! Subroutine Hypo
!-----------------------------------------------------------------------
      subroutine Hypo(stressNew,stateNew,defgradNew,
     +               stressOld,stateOld,defgradOld,dt,props,
     +               nblock,nstatev,nprops,Dissipation)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nblock, nstatev, nprops
      real*8 dt
      real*8 props(nprops),defgradOld(nblock,9),
     +  stressOld(nblock,6),
     +  stateOld(nblock,nstatev),
     +  defgradNew(nblock,9),
     +  stressNew(nblock,6), stateNew(nblock,nstatev)
!-----------------------------------------------------------------------
!     	Internal vumat variables
!-----------------------------------------------------------------------
      integer alpha,hflag,km,Txflag
      parameter(alpha=12)! Number of slip systems (12 for FCC materials)
      real*8 C11! Elastic coefficient
      real*8 C12! Elastic coefficient
      real*8 C44! Elastic coefficient
      real*8 gamma! Accumulated plastic shear strain
      real*8 PEQ! Equivalent von Mises plastic strain
      real*8 dgamma(alpha)! Shear strain increment for slip system alpha
      real*8 gamma0_dot! Referance slip rate
      real*8 bm! Instantaneous strain rate sensitivity
      real*8 tau0_c! Initial critical resolved shear stress
      real*8 theta1! Hardening parameter (Voce)
      real*8 tau1! Hardening parameter (Voce)
      real*8 theta2! Hardening parameter (Voce)
      real*8 tau2! Hardening parameter (Voce)
      real*8 h0! Hardening parameter
      real*8 tau_s! Hardening parameter
      real*8 am! Hardening parameter
      real*8 dtau_c(alpha)! Critical resolved shear stress increment
! for slip system alpha
      real*8 q(12,12)! Latent hardening matrix
      real*8 tau(alpha)! Resolved shear stress for slip system alpha
      real*8 tau_c(alpha)! Critical resolved shear stress 
! for slip system alpha
      real*8 n(alpha,3)! Slip plane normal for slip system alpha
      real*8 m(alpha,3)! Slip direction for slip system alpha
      real*8 Fold(3,3)! Old Deformation gradient F=RU
      real*8 Fnew(3,3)! New Deformation gradient F=RU
      real*8 R(3,3),RT(3,3)! Rotation tensor w.r.t. W and its transpose
      real*8 phi1, PHI, phi2! Euler angles (phi1, PHI, phi2)
      real*8 S(alpha,3,3)! Schmid tensor for slip system alpha
      integer a,i,j! Loop variables
      real*8 sigs(6)! Stress tensor components, 
! S11, S22, S33, S12, S23, S31 in global coordinate system
      real*8 sigma(6)! Corotaional stress tensor components, 
! S11, S22, S33, S12, S23, S31 w.r.t. W
      real*8 depsilon(6)! Corotaional incremental strain tensor 
! components, dE11, dE22, dE33, dE12, dE23, dE31 w.r.t. W
      real*8 depsilon_p(6)! Corotaional incremental plastic strain
! tensor components, dE11, dE22, dE33, dE12, dE23, dE31 w.r.t. W
      real*8 domega_p(3)! Corotaional incremental plastic spin tensor 
! components, dW32, dW13, dW21 w.r.t. W
      real*8 domega_e(3)! Incremental elastic spin tensor components,
! dW32, dW13, dW21 in global coordinate system
      real*8 spininc(3)! Incremental spin tensor components, 
! dW32, dW13, dW21 in global coordinate system
      real*8 epsinc(6)! Incremental incremental strain tensor 
! components, dE11, dE22, dE33, dE12, dE23, dE31 in the 
! global coordinate system
      real*8 xmat1(3,3), xmat2(3,3)! Tensors used for transformations
      real*8 Dissipation(nblock)! The change in dissipated inelastic 
! specific energy (sigma_ij*D^p_ij*dt=sum(tau(alpha)*dgamma(alpha)))
      real*8 ang(3)! Euler angles phi1, PHI, phi2
      real*8 four, three, two, one, half, zero, halfCirc
      real*8 Pi, oSqrtThree, oSqrtTwo, small, critEps
      parameter (Pi=4.d0*atan(1.d0))
      parameter(four=4.d0, three=3.d0, two=2.d0, one=1.d0,
     +          half=5d-1, zero=0.d0, oSqrtThree=1.d0/sqrt(3.d0),
     +          oSqrtTwo=1.d0/sqrt(2.d0), small=1.d-6, critEps=1.d-6,
     +          halfCirc=180.d0)! Constants
      integer nsub,k! Nuber of sub-steps and sub-step loop variable
      real*8 dti! Sub-stepping time step
#if SCMM_HYPO_DFLAG != 0
      real*8 VVF0, VVFC, VVF, q1, q2 ! Damage variables
#endif
#if SCMM_HYPO_DFLAG != 0 || SCMM_HYPO_TFLAG != 0
      integer isActive ! Is the integration point active (0=deleted, 
!                                                         1=active)
#endif
#if SCMM_HYPO_TFLAG != 0
      integer Tflag
      real*8 rho, CP, xsi, Ti, T0, Tm, Tc, mt ! Adiabatic heating params
      real*8 temperature, normTemp
#endif
      real*8 deltaDisp
!-----------------------------------------------------------------------
!     Read parameters from ABAQUS material card
!-----------------------------------------------------------------------
      C11        = props(1)! Elastic coefficient
      C12        = props(2)! Elastic coefficient
      C44        = props(3)! Elastic coefficient
      gamma0_dot = props(4)! Referance slip rate
      bm         = props(5)! Instantaneous strain rate sensitivity
      tau0_c     = props(6)! Initial critical resolved shear stress
! Texture flag (1=Euler angle from material card,
!               2=Euler angle from history card)
      Txflag     = nint(props(8))
      phi1       = props(9)*Pi/halfCirc! Euler angle phi1 in radians
      PHI        = props(10)*Pi/halfCirc! Euler angle PHI in radians
      phi2       = props(11)*Pi/halfCirc! Euler angle phi2 in radians
      hflag      = nint(props(12))! Hardening type (1=Voce,2=Kalidindi)
#if SCMM_HYPO_DFLAG != 0
      VVF0       = props(18) ! Initial damage / void volume fraction
      VVFC       = props(19) ! Critical damage / void volume fraction
      q1         = props(20) ! Damage evolution parameter
      q2         = props(21) ! Damage evolution parameter
#endif
#if SCMM_HYPO_TFLAG !=0
      Tflag      = nint(props(22)) ! Adiabatic heating (1=on, else=off)
      rho        = props(23) ! Density used in adiabatic heating
      CP         = props(24) ! Specific heat capacity
      xsi        = props(25) ! Taylor-Quinney coefficient
      Ti         = props(26) ! Initial temperature
      T0         = props(27) ! Reference temperature
      Tm         = props(28) ! Melting temperature
      Tc         = props(29) ! Critical temperature
      mt         = props(30) ! Exponent for temperature dependency
#endif
!-----------------------------------------------------------------------
!     Determine the hardening law parameters
!-----------------------------------------------------------------------
#ifdef SCMM_HYPO_VOCE_ONLY
!-----------------------------------------------------------------------
!       Voce
!-----------------------------------------------------------------------
      call unpackVoce(nprops,props,theta1,tau1,theta2,tau2,q)
#elif defined SCMM_HYPO_KALIDINDI_ONLY
!-----------------------------------------------------------------------
!       Kalidindi et al.
!-----------------------------------------------------------------------
      call unpackKalidindi(nprops,props,h0,tau_s,am,q)
#else
      if(hflag.eq.1)then
!-----------------------------------------------------------------------
!       Voce
!-----------------------------------------------------------------------
        call unpackVoce(nprops,props,theta1,tau1,theta2,tau2,q)
      elseif(hflag.eq.2)then
!-----------------------------------------------------------------------
!       Kalidindi et al.
!-----------------------------------------------------------------------
        call unpackKalidindi(nprops,props,h0,tau_s,am,q)
      else
!-----------------------------------------------------------------------
!       Error on wrong hflag
!-----------------------------------------------------------------------
#if defined SCMM_HYPO_STANDARD
         call STDB_ABQERR(-3,'Wrong Hardening model, hflag = %I',
     .                    hflag,,)
#elif defined SCMM_HYPO_EXPLICIT
         call XPLB_ABQERR(-3,'Wrong Hardening model, hflag = %I',
     .                    hflag,,)
#else
         write(*,*) 'hflag = ',hflag
         error stop 'ERROR: Wrong Hardening model'
#endif
      endif
#endif
!-----------------------------------------------------------------------
!     Slip normals and directions in local coordinate system for FCC
!-----------------------------------------------------------------------
      n(1,1:3) = (/ oSqrtThree, oSqrtThree, oSqrtThree/)
      n(2,1:3) = (/ oSqrtThree, oSqrtThree, oSqrtThree/)
      n(3,1:3) = (/ oSqrtThree, oSqrtThree, oSqrtThree/)
      n(4,1:3) = (/-oSqrtThree,-oSqrtThree, oSqrtThree/)
      n(5,1:3) = (/-oSqrtThree,-oSqrtThree, oSqrtThree/)
      n(6,1:3) = (/-oSqrtThree,-oSqrtThree, oSqrtThree/)
      n(7,1:3) = (/-oSqrtThree, oSqrtThree, oSqrtThree/)
      n(8,1:3) = (/-oSqrtThree, oSqrtThree, oSqrtThree/)
      n(9,1:3) = (/-oSqrtThree, oSqrtThree, oSqrtThree/)
      n(10,1:3)= (/ oSqrtThree,-oSqrtThree, oSqrtThree/)
      n(11,1:3)= (/ oSqrtThree,-oSqrtThree, oSqrtThree/)
      n(12,1:3)= (/ oSqrtThree,-oSqrtThree, oSqrtThree/)
!-----------------------------------------------------------------------
      m(1,1:3) = (/-oSqrtTwo, zero    , oSqrtTwo/)
      m(2,1:3) = (/-oSqrtTwo, oSqrtTwo, zero    /)
      m(3,1:3) = (/ zero    ,-oSqrtTwo, oSqrtTwo/)
      m(4,1:3) = (/ zero    , oSqrtTwo, oSqrtTwo/)
      m(5,1:3) = (/-oSqrtTwo, oSqrtTwo, zero    /)
      m(6,1:3) = (/ oSqrtTwo, zero    , oSqrtTwo/)
      m(7,1:3) = (/ oSqrtTwo, zero    , oSqrtTwo/)
      m(8,1:3) = (/ oSqrtTwo, oSqrtTwo, zero    /)
      m(9,1:3) = (/ zero    ,-oSqrtTwo, oSqrtTwo/)
      m(10,1:3)= (/ zero    , oSqrtTwo, oSqrtTwo/)
      m(11,1:3)= (/ oSqrtTwo, oSqrtTwo, zero    /)
      m(12,1:3)= (/-oSqrtTwo, zero    , oSqrtTwo/)
!-----------------------------------------------------------------------
      do j=1,3
         do i=1,3
            do a=1,alpha
               S(a,i,j)=m(a,i)*n(a,j)
            enddo
         enddo
      enddo
!-----------------------------------------------------------------------
!     Time greater than zero
!-----------------------------------------------------------------------
      if(stateold(1,13).lt.small)then ! First step
!-----------------------------------------------------------------------
!       Initializing the rotation tensor
!-----------------------------------------------------------------------
        if (Txflag.eq.2)then
!-----------------------------------------------------------------------
!         Load orientations from initial conditions
!-----------------------------------------------------------------------
          do km=1,nblock
            phi1 = STATEOLD(km,1)*Pi/halfCirc
            PHI  = STATEOLD(km,2)*Pi/halfCirc
            phi2 = STATEOLD(km,3)*Pi/halfCirc
!-----------------------------------------------------------------------
            R(1,1) =  cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(PHI)
            R(1,2) = -cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(PHI)
            R(1,3) =  sin(phi1)*sin(PHI)
            R(2,1) =  sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(PHI)
            R(2,2) = -sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(PHI)
            R(2,3) = -cos(phi1)*sin(PHI)
            R(3,1) =  sin(phi2)*sin(PHI)
            R(3,2) =  cos(phi2)*sin(PHI)
            R(3,3) =  cos(PHI)
!-----------------------------------------------------------------------
            a = 4
            do j=1,3
              do i=1,3
                STATEOLD(km,a) = R(i,j)
                a              = a+1
              enddo
            enddo
          enddo
        else
!-----------------------------------------------------------------------
!         Load orientation from material properties
!-----------------------------------------------------------------------
          R(1,1) =  cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(PHI)
          R(1,2) = -cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(PHI)
          R(1,3) =  sin(phi1)*sin(PHI)
          R(2,1) =  sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(PHI)
          R(2,2) = -sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(PHI)
          R(2,3) = -cos(phi1)*sin(PHI)
          R(3,1) =  sin(phi2)*sin(PHI)
          R(3,2) =  cos(phi2)*sin(PHI)
          R(3,3) =  cos(PHI)
!-----------------------------------------------------------------------
          a = 4
          do j=1,3
            do i=1,3
              do km=1,nblock
                STATEOLD(km,a) = R(i,j)
              enddo
              a              = a+1
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!       Initializing the other state variables
!-----------------------------------------------------------------------
        do km=1,nblock
          STATEOLD(km,13:24) = tau0_c
          STATEOLD(km,25)    = zero
          STATEOLD(km,27)    = zero
#if SCMM_HYPO_DFLAG != 0
          STATEOLD(km,29)    = VVF0
#endif
#if SCMM_HYPO_DFLAG != 0 || SCMM_HYPO_TFLAG != 0
          STATEOLD(km,30)    = one
#endif
        enddo
#if SCMM_HYPO_TFLAG == 1
        do km=1,nblock
          STATEOLD(km,31)  = Ti
        enddo
#elif SCMM_HYPO_TFLAG != 0
        if(Tflag.eq.1)then
          do km=1,nblock
            STATEOLD(km,31)  = Ti
          enddo
        endif
#endif
      endif
!-----------------------------------------------------------------------
!     Loop over nblock integration points
!-----------------------------------------------------------------------
      do km = 1, nblock
!-----------------------------------------------------------------------
!       Defining state variables from last increment
!-----------------------------------------------------------------------
        a = 4
        do j=1,3
          do i=1,3
            R(i,j) = STATEOLD(km,a)
            a      = a+1
          enddo
        enddo
        tau_c = STATEOLD(km,13:24)
        gamma = STATEOLD(km,25)
        PEQ   = STATEOLD(km,27)
#if SCMM_HYPO_DFLAG != 0
        VVF   = STATEOLD(km,29)
#endif
#if SCMM_HYPO_DFLAG != 0 || SCMM_HYPO_TFLAG != 0
        isActive = nint(STATEOLD(km,30))
#endif
#if SCMM_HYPO_TFLAG == 1
        temperature = STATEOLD(km,31)
#elif SCMM_HYPO_TFLAG != 0
        if(Tflag.eq.1)then
          temperature = STATEOLD(km,31)
        endif
#endif
!-----------------------------------------------------------------------
!       Check if integration point is active
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG != 0 || SCMM_HYPO_TFLAG != 0
#ifdef SCMM_HYPO_EXPLICIT
        if(isActive.eq.0)then
            stressNew(km,1:6) = zero
            STATENEW(km,1:nstatev) = STATEOLD(km,1:nstatev)
            Dissipation(km) = zero
            cycle ! Continue to next loop cycle
        endif
#endif
#endif
!-----------------------------------------------------------------------
!       Co-rotating the stress tensor
!-----------------------------------------------------------------------
        sigs(1) = stressOld(km,1)
        sigs(2) = stressOld(km,2)
        sigs(3) = stressOld(km,3)
        sigs(4) = stressOld(km,4)
        sigs(5) = stressOld(km,5)
        sigs(6) = stressOld(km,6)
!-----------------------------------------------------------------------
!       Calculating the transpose of the rotation tensor
!-----------------------------------------------------------------------
        call mtransp(R,RT)
!-----------------------------------------------------------------------
!       Stress components, sigma_hat=R**T sigma R
!-----------------------------------------------------------------------
        call vec2mat(sigs,xmat1)
        call transform(xmat1,RT,R,xmat2)
        call mat2vec(xmat2,sigma)
!-----------------------------------------------------------------------
!       Calculating the effective stress sigma_eff=sigma/(1-VVF)
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG != 0
        sigma = sigma/(one-VVF)
#endif
!-----------------------------------------------------------------------
!       Calculating the strain and spin increments from
!       the deformation gradient in the global coordinate system
!-----------------------------------------------------------------------
!       Old deformation gradient, F
!-----------------------------------------------------------------------
        Fold(1,1) = defgradOld(km,1)
        Fold(2,2) = defgradOld(km,2)
        Fold(3,3) = defgradOld(km,3)
        Fold(1,2) = defgradOld(km,4)
        Fold(2,3) = defgradOld(km,5)
        Fold(3,1) = defgradOld(km,6)
        Fold(2,1) = defgradOld(km,7)
        Fold(3,2) = defgradOld(km,8)
        Fold(1,3) = defgradOld(km,9)
!-----------------------------------------------------------------------
!       New deformation gradient, F
!-----------------------------------------------------------------------
        Fnew(1,1) = defgradNew(km,1)
        Fnew(2,2) = defgradNew(km,2)
        Fnew(3,3) = defgradNew(km,3)
        Fnew(1,2) = defgradNew(km,4)
        Fnew(2,3) = defgradNew(km,5)
        Fnew(3,1) = defgradNew(km,6)
        Fnew(2,1) = defgradNew(km,7)
        Fnew(3,2) = defgradNew(km,8)
        Fnew(1,3) = defgradNew(km,9)
        call sinc(Fold,Fnew,dt,epsinc,spininc)
!-----------------------------------------------------------------------
!       Begin the sub-stepping
!-----------------------------------------------------------------------
        nsub = ceiling(sqrt(epsinc(1)**two+epsinc(2)**two+
     +                      epsinc(3)**two+two*epsinc(4)**two+
     +                      two*epsinc(5)**two+
     +                      two*epsinc(6)**two)/(critEps))
!-----------------------------------------------------------------------
        epsinc  = epsinc/nsub
        spininc = spininc/nsub
        dti     = dt/nsub
        Dissipation(km) = zero
!-----------------------------------------------------------------------
        do k=1,nsub
!-----------------------------------------------------------------------
!       Strain increments, depsilon_hat=R**T depsilon R
!-----------------------------------------------------------------------
          call vec2mat(epsinc,xmat1)
          call transform(xmat1,RT,R,xmat2)
          call mat2vec(xmat2,depsilon)
!-----------------------------------------------------------------------
!         Calculating resolved shear stresses at n and estimate the 
!         shear strain increments at n+1
!-----------------------------------------------------------------------
          depsilon_p = zero
          domega_p   = zero
          dgamma     = zero
          dtau_c     = zero
          deltaDisp  = zero
!-----------------------------------------------------------------------
#if SCMM_HYPO_TFLAG == 1
          normTemp = ((temperature-T0)/(Tm-T0))**mt
#elif SCMM_HYPO_TFLAG != 0
          if(Tflag.eq.1)then
            normTemp = ((temperature-T0)/(Tm-T0))**mt
          else
            normTemp = zero
          endif
#endif
!-----------------------------------------------------------------------
          do a=1,alpha
            tau(a) = sigma(1)*S(a,1,1)+sigma(2)*S(a,2,2)+
     +               sigma(3)*S(a,3,3)+sigma(4)*(S(a,1,2)+S(a,2,1))+
     +               sigma(5)*(S(a,2,3)+S(a,3,2))+
     +               sigma(6)*(S(a,3,1)+S(a,1,3))
#if SCMM_HYPO_TFLAG == 0
            dgamma(a) = dti*gamma0_dot*
     +                (abs(tau(a)/tau_c(a)))**(one/bm)*sign(one,tau(a))
#else
            dgamma(a) = dti*gamma0_dot*
     +                (abs(tau(a)/(tau_c(a)*(one-normTemp))))
     +                **(one/bm)*sign(one,tau(a))
#endif
!-----------------------------------------------------------------------
!       Calculating corotated incremental plastic strain and spin
!-----------------------------------------------------------------------
            depsilon_p(1) = depsilon_p(1)+dgamma(a)*S(a,1,1)
            depsilon_p(2) = depsilon_p(2)+dgamma(a)*S(a,2,2)
            depsilon_p(3) = depsilon_p(3)+dgamma(a)*S(a,3,3)
            depsilon_p(4) = depsilon_p(4)+
     +                      half*dgamma(a)*(S(a,1,2)+S(a,2,1))
            depsilon_p(5) = depsilon_p(5)+
     +                      half*dgamma(a)*(S(a,2,3)+S(a,3,2))
            depsilon_p(6) = depsilon_p(6)+
     +                      half*dgamma(a)*(S(a,3,1)+S(a,1,3))
!-----------------------------------------------------------------------
            domega_p(1) = domega_p(1)+
     +                    half*dgamma(a)*(S(a,3,2)-S(a,2,3))
            domega_p(2) = domega_p(2)+
     +                    half*dgamma(a)*(S(a,1,3)-S(a,3,1))
            domega_p(3) = domega_p(3)+
     +                    half*dgamma(a)*(S(a,2,1)-S(a,1,2))
!-----------------------------------------------------------------------
!       Approximating the dissipated energy by using tau at n
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG != 0
            deltaDisp = deltaDisp+(one-VVF)*tau(a)*dgamma(a)
#else
            deltaDisp = deltaDisp+tau(a)*dgamma(a)
#endif
          enddo
          Dissipation(km) = Dissipation(km)+deltaDisp
!-----------------------------------------------------------------------
!       Updating temperature due to adiabatic heating
!-----------------------------------------------------------------------
#if SCMM_HYPO_TFLAG == 1
          temperature = temperature+xsi*deltaDisp/(rho*CP)
#elif SCMM_HYPO_TFLAG != 0
          if(Tflag.eq.1)then
            temperature = temperature+xsi*deltaDisp/(rho*CP)
          endif
#endif
!-----------------------------------------------------------------------
!       Updating corotated stress tensor
!-----------------------------------------------------------------------
          sigma(1) = sigma(1)+C11*(depsilon(1)-depsilon_p(1))
     +                       +C12*(depsilon(2)-depsilon_p(2))
     +                       +C12*(depsilon(3)-depsilon_p(3))
          sigma(2) = sigma(2)+C12*(depsilon(1)-depsilon_p(1))
     +                       +C11*(depsilon(2)-depsilon_p(2))
     +                       +C12*(depsilon(3)-depsilon_p(3))
          sigma(3) = sigma(3)+C12*(depsilon(1)-depsilon_p(1))
     +                       +C12*(depsilon(2)-depsilon_p(2))
     +                       +C11*(depsilon(3)-depsilon_p(3))
          sigma(4) = sigma(4)+two*C44*(depsilon(4)-depsilon_p(4))
          sigma(5) = sigma(5)+two*C44*(depsilon(5)-depsilon_p(5))
          sigma(6) = sigma(6)+two*C44*(depsilon(6)-depsilon_p(6))
!-----------------------------------------------------------------------
!       Updating critical resolved shear stresses
!-----------------------------------------------------------------------
#ifdef SCMM_HYPO_VOCE_ONLY
          call Voce(q,theta1,tau1,theta2,
     +              tau2,dgamma,gamma,tau_c)
#elif defined SCMM_HYPO_KALIDINDI_ONLY
          call Kalidindi(q,h0,tau_s,am,dgamma,tau_c)
#else
          if(hflag.eq.1)then
            call Voce(q,theta1,tau1,theta2,
     +                tau2,dgamma,gamma,tau_c)
! hflag=2 (have already checked if hflag is not equal to 1 or 2)
          else 
            call Kalidindi(q,h0,tau_s,am,dgamma,tau_c)
          endif
#endif
!-----------------------------------------------------------------------
!       Updating accumulated plastic shear strain
!-----------------------------------------------------------------------
          gamma =
     +         gamma+abs(dgamma(1))+abs(dgamma(2))+abs(dgamma(3))+
     +         abs(dgamma(4))+abs(dgamma(5))+abs(dgamma(6))+
     +         abs(dgamma(7))+abs(dgamma(8))+abs(dgamma(9))+
     +         abs(dgamma(10))+abs(dgamma(11))+abs(dgamma(12))
!-----------------------------------------------------------------------
!       Equivalent von mises plastic strain
!-----------------------------------------------------------------------
          PEQ = PEQ+sqrt(two*(depsilon_p(1)**two+
     +                    depsilon_p(2)**two+depsilon_p(3)**two+
     +                    two*depsilon_p(4)**two+
     +                    two*depsilon_p(5)**two+
     +                    two*depsilon_p(6)**two)/three)
!-----------------------------------------------------------------------
!       Calculating incremental elastic rotation in the 
!       global coordinate system
!-----------------------------------------------------------------------
          xmat1(1,1) = zero
          xmat1(1,2) = -domega_p(3)
          xmat1(1,3) = domega_p(2)
          xmat1(2,1) = domega_p(3)
          xmat1(2,2) = zero
          xmat1(2,3) = -domega_p(1)
          xmat1(3,1) = -domega_p(2)
          xmat1(3,2) = domega_p(1)
          xmat1(3,3) = zero
!-----------------------------------------------------------------------
          call transform(xmat1,R,RT,xmat2)
!-----------------------------------------------------------------------
          domega_e(1) = spininc(1)-xmat2(3,2)
          domega_e(2) = spininc(2)-xmat2(1,3)
          domega_e(3) = spininc(3)-xmat2(2,1)
!-----------------------------------------------------------------------
!       Updating the rotation tensor
!-----------------------------------------------------------------------
          call updateR(domega_e,R)
          call mtransp(R,RT)
!-----------------------------------------------------------------------
!       Updating damage
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG != 0
          call UpdateDamage(VVF,sigma,dgamma,q1,q2)
#endif
!-----------------------------------------------------------------------
!       End sub-stepping
!-----------------------------------------------------------------------
        enddo! End sub-stepping
!-----------------------------------------------------------------------
!       Check for failure
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG != 0
        if((VVF.ge.VVFC).or.(VVF.ge.one))then
          VVF = min(VVFC,one)
          isActive = 0
        endif
#endif
#if SCMM_HYPO_TFLAG == 1
        if((temperature.ge.Tc).or.(temperature.ge.Tm))then
          temperature = min(temperature,Tm)
          isActive = 0
        endif
#elif SCMM_HYPO_TFLAG != 0
        if(Tflag.eq.1)then
          if((temperature.ge.Tc).or.(temperature.ge.Tm))then
            temperature = min(temperature,Tm)
            isActive = 0
          endif
        endif
#endif
!-----------------------------------------------------------------------
!       Calculating the Cauchy stress tensor from the effective stress
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG != 0
        sigma = sigma*(one-VVF)
#endif
!-----------------------------------------------------------------------
!       Transform the stress tensor back to the global coordinate system
!-----------------------------------------------------------------------
        call vec2mat(sigma,xmat1)
        call transform(xmat1,R,RT,xmat2)
        call mat2vec(xmat2,sigs)
        stressNew(km,1) = sigs(1)
        stressNew(km,2) = sigs(2)
        stressNew(km,3) = sigs(3)
        stressNew(km,4) = sigs(4)
        stressNew(km,5) = sigs(5)
        stressNew(km,6) = sigs(6)
!-----------------------------------------------------------------------
!       Updating output variables
!-----------------------------------------------------------------------
        a = 4
        do j=1,3
          do i=1,3
            STATENEW(km,a) = R(i,j)! Rotation tensor
            a              = a+1
          enddo
        enddo
        ! Critical resolved shear stresses/ Slip resistances
        STATENEW(km,13:24) = tau_c
        STATENEW(km,25)    = gamma! Accumulated plastic strain
        ! Equivalent von Mises stress
        STATENEW(km,26)    = sqrt(half*((sigma(1)-sigma(2))**two
     +                                 +(sigma(2)-sigma(3))**two
     +                                 +(sigma(3)-sigma(1))**two)
     +                      +three*sigma(4)**two+three*sigma(5)**two
     +                      +three*sigma(6)**two)
        STATENEW(km,27) = PEQ! Equivalent von mises plastic strain
        STATENEW(km,28) = nsub! Number of sub steps
#if SCMM_HYPO_DFLAG != 0
        STATENEW(km,29) = VVF ! Damage / void volume fraction
#endif
! Is the element active or should it be deleted (Abaqus status variable)
#if SCMM_HYPO_DFLAG != 0 || SCMM_HYPO_TFLAG != 0
        STATENEW(km,30) = isActive
#endif
#if SCMM_HYPO_TFLAG == 1
        STATENEW(km,31) = temperature
#elif SCMM_HYPO_TFLAG != 0
        if(Tflag.eq.1)then
          STATENEW(km,31) = temperature
        endif
#endif
!-----------------------------------------------------------------------
        call euler(R,ang)
!-----------------------------------------------------------------------
        STATENEW(km,1:3) = ang! Euler angles phi1, PHI, phi2
!-----------------------------------------------------------------------
!       end loops
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      return
      end subroutine Hypo
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------