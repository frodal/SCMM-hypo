!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Subroutine CCCP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Subroutines should be inlined by the compiler
!-----------------------------------------------------------------------
!DIR$ ATTRIBUTES FORCEINLINE :: CCCP
!-----------------------------------------------------------------------
! Preprocessor definitions
!-----------------------------------------------------------------------
#ifndef SCMM_HYPO_CCCP
#define SCMM_HYPO_CCCP
!-----------------------------------------------------------------------
! Subroutine CCCP
!-----------------------------------------------------------------------
      subroutine CCCP(stressNew,stateNew,defgradNew,
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
!       Internal vumat variables
!-----------------------------------------------------------------------
      integer alpha,hflag,km,Txflag,iter,maxIter
      parameter(alpha=12)! Number of slip systems (12 for FCC materials)
      real*8 C11! Elastic coefficient
      real*8 C12! Elastic coefficient
      real*8 C44! Elastic coefficient
      real*8 gamma,gamma_old! Accumulated plastic shear strain
      real*8 PEQ,PEQ_old! Equivalent von Mises plastic strain
      real*8 dgamma(alpha)! Shear strain increment for slip system alpha
      real*8 tau0_c! Initial critical resolved shear stress
      real*8 theta1! Hardening parameter (Voce)
      real*8 tau1! Hardening parameter (Voce)
      real*8 theta2! Hardening parameter (Voce)
      real*8 tau2! Hardening parameter (Voce)
      real*8 h0! Hardening parameter
      real*8 tau_s! Hardening parameter
      real*8 am! Hardening parameter
      real*8 dtau_c(alpha)! Critical resolved shear stress increment for slip system alpha
      real*8 q(12,12)! Latent hardening matrix
      real*8 tau(alpha)! Resolved shear stress for slip system alpha
#if SCMM_HYPO_DFLAG == 2
      real*8 tau_eff(alpha)
#endif
      real*8 tau_c(alpha)! Critical resolved shear stress for slip system alpha
      real*8 n(alpha,3)! Slip plane normal for slip system alpha
      real*8 m(alpha,3)! Slip direction for slip system alpha
      real*8 Fold(3,3)! Old Deformation gradient F=RU
      real*8 Fnew(3,3)! New Deformation gradient F=RU
      real*8 R(3,3),RT(3,3)! Rotation tensor w.r.t. W and its transpose
      real*8 phi1, PHI, phi2! Euler angles (phi1, PHI, phi2)
      real*8 S(alpha,3,3)! Schmid tensor for slip system alpha
      integer a,i,j! Loop variables
      real*8 sigs(6)! Stress tensor components, S11, S22, S33, S12, S23, S31 in global coordinate system
      real*8 sigma(6)! Corotaional stress tensor components, S11, S22, S33, S12, S23, S31 w.r.t. W
      real*8 sig_tr(6)! Trial corotaional stress tensor components, S11, S22, S33, S12, S23, S31 w.r.t. W
      real*8 depsilon(6)! Corotaional incremental strain tensor components, dE11, dE22, dE33, dE12, dE23, dE31 w.r.t. W
      real*8 depsilon_p(6)! Corotaional incremental plastic strain tensor components, dE11, dE22, dE33, dE12, dE23, dE31 w.r.t. W
      real*8 domega_p(3)! Corotaional incremental plastic spin tensor components, dW32, dW13, dW21 w.r.t. W
      real*8 domega_e(3)! Incremental elastic spin tensor components, dW32, dW13, dW21 in global coordinate system
      real*8 spininc(3)! Incremental spin tensor components, dW32, dW13, dW21 in global coordinate system
      real*8 epsinc(6)! Incremental incremental strain tensor components, dE11, dE22, dE33, dE12, dE23, dE31 in global coordinate system
      real*8 xmat1(3,3), xmat2(3,3)! Tensors used for transformations
      real*8 Dissipation(nblock)! The change in dissipated inelastic specific energy (sigma_ij*D^p_ij*dt=sum(tau(alpha)*dgamma(alpha)))
      real*8 ang(3)! Euler angles phi1, PHI, phi2
      real*8 four, three, two, one, half, zero, deg2rad
      real*8 oSqrtThree, oSqrtTwo, small, critEps
      real*8 rhoParameter,mParameter,f,tol
      real*8 dfdtau(alpha),dfdtau_c(alpha),dfdsigma(6)
#if SCMM_HYPO_DFLAG == 2
      real*8 dfdVVF, dfdsigmaskew(3)
#endif
      real*8 hMatrix(alpha,alpha),dlambda,ddgamma(alpha)
      parameter(four=4.d0, three=3.d0, two=2.d0, one=1.d0,
     +          half=5d-1, zero=0.d0,maxIter=1000,tol=1.d-8,
     +          oSqrtThree=1.d0/sqrt(3.d0),
     +          deg2rad=4.d0*atan(1.d0)/180.d0,
     +          oSqrtTwo=1.d0/sqrt(2.d0), small=1.d-6, critEps=1.d-6)! Constants
      integer nsub,k! Nuber of sub-steps and sub-step loop variable
      real*8 dti! Sub-stepping time step
#if SCMM_HYPO_DFLAG == 1 || SCMM_HYPO_DFLAG == 2
      real*8 VVF0, VVFC, VVF, q1, q2 ! Damage variables
      integer isActive ! Is the integration point active (0=deleted, 
!                                                         1=active)
#endif
#if SCMM_HYPO_DFLAG == 2
      real*8 aParam
#endif
!-----------------------------------------------------------------------
!     Read parameters from ABAQUS material card
!-----------------------------------------------------------------------
      C11          = props(1)! Elastic coefficient
      C12          = props(2)! Elastic coefficient
      C44          = props(3)! Elastic coefficient
      mParameter   = props(4)! 
      rhoParameter = props(5)! 
      tau0_c       = props(6)! Initial critical resolved shear stress
! Texture flag (1=Euler angle from material card,
!               2=Euler angle from history card)
      Txflag       = nint(props(8))
      phi1         = props(9)*deg2rad! Euler angle phi1 in radians
      PHI          = props(10)*deg2rad! Euler angle PHI in radians
      phi2         = props(11)*deg2rad! Euler angle phi2 in radians
      hflag        = nint(props(12))! Hardening type (1=Voce,2=Kalidindi)
#if SCMM_HYPO_DFLAG == 1 || SCMM_HYPO_DFLAG == 2
      VVF0         = props(18) ! Initial damage / void volume fraction
      VVFC         = props(19) ! Critical damage / void volume fraction
      q1           = props(20) ! Damage evolution parameter
      q2           = props(21) ! Damage evolution parameter
#endif
#if SCMM_HYPO_DFLAG == 2
      aParam       = props(22) ! Damage evolution parameter
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
        if (Txflag.eq.3)then
            call RandomTexture(stateOld,nblock,nstatev)
        elseif (Txflag.eq.2)then
!-----------------------------------------------------------------------
!         Load orientations from initial conditions
!-----------------------------------------------------------------------
          do km=1,nblock
            phi1 = STATEOLD(km,1)*deg2rad
            PHI  = STATEOLD(km,2)*deg2rad
            phi2 = STATEOLD(km,3)*deg2rad
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
#if SCMM_HYPO_DFLAG == 1 || SCMM_HYPO_DFLAG == 2
          STATEOLD(km,29)    = VVF0
          STATEOLD(km,30)    = one
#endif
        enddo
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
#if SCMM_HYPO_DFLAG == 1 || SCMM_HYPO_DFLAG == 2
        VVF   = STATEOLD(km,29)
        isActive = nint(STATEOLD(km,30))
!-----------------------------------------------------------------------
!       Check if integration point is active
!-----------------------------------------------------------------------
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
#if SCMM_HYPO_DFLAG == 1
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
        nsub = ceiling(sqrt(epsinc(1)**2+epsinc(2)**2+
     +                      epsinc(3)**2+two*epsinc(4)**2+
     +                      two*epsinc(5)**2+
     +                      two*epsinc(6)**2)/(critEps))
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
!       Initialize plastic variables for this sub increment
!-----------------------------------------------------------------------
          depsilon_p = zero
          domega_p   = zero
          dgamma     = zero
          dtau_c     = zero
          iter       = 0
          gamma_old  = gamma
          PEQ_old    = PEQ
!-----------------------------------------------------------------------
!       Elastic predictor (Trial stress)
!-----------------------------------------------------------------------
          sig_tr(1) = sigma(1)+C11*(depsilon(1))
     +                        +C12*(depsilon(2))
     +                        +C12*(depsilon(3))
          sig_tr(2) = sigma(2)+C12*(depsilon(1))
     +                        +C11*(depsilon(2))
     +                        +C12*(depsilon(3))
          sig_tr(3) = sigma(3)+C12*(depsilon(1))
     +                        +C12*(depsilon(2))
     +                        +C11*(depsilon(3))
          sig_tr(4) = sigma(4)+two*C44*(depsilon(4))
          sig_tr(5) = sigma(5)+two*C44*(depsilon(5))
          sig_tr(6) = sigma(6)+two*C44*(depsilon(6))
          sigma = sig_tr
!-----------------------------------------------------------------------
!       Calculating resolved shear stress for the trial state
!-----------------------------------------------------------------------
          do a=1,alpha
            tau(a) = sigma(1)*S(a,1,1)+sigma(2)*S(a,2,2)+
     +               sigma(3)*S(a,3,3)+sigma(4)*(S(a,1,2)+S(a,2,1))+
     +               sigma(5)*(S(a,2,3)+S(a,3,2))+
     +               sigma(6)*(S(a,3,1)+S(a,1,3))
          enddo
#if SCMM_HYPO_DFLAG == 2
          call calcTauEff(tau,sigma,VVF,aParam,q1,q2,tau_eff)
#endif
!-----------------------------------------------------------------------
!       Calculate the yield function based on the trial state
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 2
          call yieldfunction(tau_eff,tau_c,
     .                       rhoParameter,mParameter,f)
#else
          call yieldfunction(tau,tau_c,rhoParameter,mParameter,f)
#endif
!-----------------------------------------------------------------------
!       Check yield criterion
!-----------------------------------------------------------------------
          if (f.gt.zero)then
!-----------------------------------------------------------------------
!       Return mapping (Cutting plane)
!-----------------------------------------------------------------------
            do while ((abs(f).gt.tol).and.(iter.lt.maxIter))
              iter = iter+1
!-----------------------------------------------------------------------
!       Calculating gradients to the yield function
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 2
              call yieldgradient(tau,tau_eff,sigma,tau_c,VVF,
     .                           rhoParameter,mParameter,aParam,
     .                           q1,q2,S,dfdtau,dfdtau_c,dfdsigma,
     .                           dfdsigmaskew,dfdVVF)
#else
              call yieldgradient(tau,tau_c,rhoParameter,mParameter,
     .                           S,dfdtau,dfdtau_c,dfdsigma)
#endif
!-----------------------------------------------------------------------
!       Calculating the work-hardening rate matrix
!-----------------------------------------------------------------------
#ifdef SCMM_HYPO_VOCE_ONLY
              call VoceMatrix(q,theta1,tau1,theta2,tau2,gamma,hMatrix)
#elif defined SCMM_HYPO_KALIDINDI_ONLY
              call KalidindiMatrix(q,h0,tau_s,am,tau_c,hMatrix)
#else
              if(hflag.eq.1)then
                call VoceMatrix(q,theta1,tau1,theta2,
     .                          tau2,gamma,hMatrix)
              else ! hflag=2 (have already checked if hflag is not equal to 1 or 2)
                call KalidindiMatrix(q,h0,tau_s,am,
     .                               tau_c,hMatrix)
              endif
#endif
!-----------------------------------------------------------------------
!       Calculate increment in plastic parameter deltaLambda
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 2
              call RMAP(f,dfdtau,dfdtau_c,dfdsigma,VVF,dfdVVF,
     .                  sigma,C11,C12,C44,hMatrix,dlambda)
#else
              call RMAP(f,dfdtau,dfdtau_c,
     .                  dfdsigma,C11,C12,C44,hMatrix,dlambda)
#endif
!-----------------------------------------------------------------------
!       Update plastic slip dgamma(alpha)
!-----------------------------------------------------------------------
              ddgamma = dlambda*dfdtau
              dgamma = dgamma+ddgamma
!-----------------------------------------------------------------------
!       Updating accumulated plastic shear strain
!-----------------------------------------------------------------------
              gamma = gamma_old+
     +                abs(dgamma(1))+abs(dgamma(2))+abs(dgamma(3))+
     +                abs(dgamma(4))+abs(dgamma(5))+abs(dgamma(6))+
     +                abs(dgamma(7))+abs(dgamma(8))+abs(dgamma(9))+
     +                abs(dgamma(10))+abs(dgamma(11))+abs(dgamma(12))
!-----------------------------------------------------------------------
!       Update plastic strain increment and plastic spin increment
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 2
              depsilon_p(1) = depsilon_p(1)+
     +                        dlambda*dfdsigma(1)*(one-VVF)
              depsilon_p(2) = depsilon_p(2)+
     +                        dlambda*dfdsigma(2)*(one-VVF)
              depsilon_p(3) = depsilon_p(3)+
     +                        dlambda*dfdsigma(3)*(one-VVF)
              depsilon_p(4) = depsilon_p(4)+
     +                        dlambda*dfdsigma(4)*(one-VVF)
              depsilon_p(5) = depsilon_p(5)+
     +                        dlambda*dfdsigma(5)*(one-VVF)
              depsilon_p(6) = depsilon_p(6)+
     +                        dlambda*dfdsigma(6)*(one-VVF)
!
              domega_p(1) = domega_p(1)+
     +                      dlambda*dfdsigmaskew(1)*(one-VVF)
              domega_p(2) = domega_p(2)+
     +                      dlambda*dfdsigmaskew(2)*(one-VVF)
              domega_p(3) = domega_p(3)+
     +                      dlambda*dfdsigmaskew(3)*(one-VVF)
#else
              do a=1,alpha
                depsilon_p(1) = depsilon_p(1)+ddgamma(a)*S(a,1,1)
                depsilon_p(2) = depsilon_p(2)+ddgamma(a)*S(a,2,2)
                depsilon_p(3) = depsilon_p(3)+ddgamma(a)*S(a,3,3)
                depsilon_p(4) = depsilon_p(4)+
     +                          half*ddgamma(a)*(S(a,1,2)+S(a,2,1))
                depsilon_p(5) = depsilon_p(5)+
     +                          half*ddgamma(a)*(S(a,2,3)+S(a,3,2))
                depsilon_p(6) = depsilon_p(6)+
     +                          half*ddgamma(a)*(S(a,3,1)+S(a,1,3))
!
                domega_p(1) = domega_p(1)+
     +                        half*ddgamma(a)*(S(a,3,2)-S(a,2,3))
                domega_p(2) = domega_p(2)+
     +                        half*ddgamma(a)*(S(a,1,3)-S(a,3,1))
                domega_p(3) = domega_p(3)+
     +                        half*ddgamma(a)*(S(a,2,1)-S(a,1,2))
              enddo
#endif
!-----------------------------------------------------------------------
!       Equivalent von mises plastic strain
!-----------------------------------------------------------------------
              PEQ = PEQ_old+sqrt(two*(depsilon_p(1)**2+
     +                      depsilon_p(2)**2+depsilon_p(3)**2+
     +                      two*depsilon_p(4)**2+
     +                      two*depsilon_p(5)**2+
     +                      two*depsilon_p(6)**2)/three)
!-----------------------------------------------------------------------
!       Updating damage
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 2
              call UpdateDamageHan(VVF,dfdsigma,dlambda)
              if((VVF.ge.VVFC).or.(VVF.ge.one))then
                VVF = min(VVFC,one)
                isActive = 0
              endif
#endif
!-----------------------------------------------------------------------
!       Updating corotated stress tensor and resolved shear stress
!-----------------------------------------------------------------------
              sigma(1) = sig_tr(1)+C11*(-depsilon_p(1))
     +                            +C12*(-depsilon_p(2))
     +                            +C12*(-depsilon_p(3))
              sigma(2) = sig_tr(2)+C12*(-depsilon_p(1))
     +                            +C11*(-depsilon_p(2))
     +                            +C12*(-depsilon_p(3))
              sigma(3) = sig_tr(3)+C12*(-depsilon_p(1))
     +                            +C12*(-depsilon_p(2))
     +                            +C11*(-depsilon_p(3))
              sigma(4) = sig_tr(4)+two*C44*(-depsilon_p(4))
              sigma(5) = sig_tr(5)+two*C44*(-depsilon_p(5))
              sigma(6) = sig_tr(6)+two*C44*(-depsilon_p(6))
!-----------------------------------------------------------------------
              do a=1,alpha
                tau(a) = sigma(1)*S(a,1,1)+sigma(2)*S(a,2,2)+
     +                   sigma(3)*S(a,3,3)+sigma(4)*(S(a,1,2)+S(a,2,1))+
     +                   sigma(5)*(S(a,2,3)+S(a,3,2))+
     +                   sigma(6)*(S(a,3,1)+S(a,1,3))
              enddo
#if SCMM_HYPO_DFLAG == 2
              call calcTauEff(tau,sigma,VVF,aParam,q1,q2,tau_eff)
#endif
!-----------------------------------------------------------------------
!       Updating critical resolved shear stresses
!-----------------------------------------------------------------------
#ifdef SCMM_HYPO_VOCE_ONLY
              call VoceCCCP(q,theta1,tau1,theta2,
     +                      tau2,dfdtau,dlambda,gamma,tau_c)
#elif defined SCMM_HYPO_KALIDINDI_ONLY
              call KalidindiCCCP(q,h0,tau_s,am,
     +                           dfdtau,dlambda,tau_c)
#else
              if(hflag.eq.1)then
                call VoceCCCP(q,theta1,tau1,theta2,
     +                        tau2,dfdtau,dlambda,gamma,tau_c)
              else ! hflag=2 (have already checked if hflag is not equal to 1 or 2)
                call KalidindiCCCP(q,h0,tau_s,am,
     +                             dfdtau,dlambda,tau_c)
              endif
#endif
!-----------------------------------------------------------------------
!       Update the yield function
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 2
              call yieldfunction(tau_eff,tau_c,
     .                             rhoParameter,mParameter,f)
#else
              call yieldfunction(tau,tau_c,rhoParameter,mParameter,f)
#endif
            enddo
            if ((iter.ge.maxIter).and.(abs(f).gt.tol))then
#if defined SCMM_HYPO_STANDARD
         call STDB_ABQERR(-3,'Maximum number of RMAP iterations'//
     . ' reached. Maximum number of iterations: %I, abs(f) = %R',
     .                    maxIter,abs(f),)
#elif defined SCMM_HYPO_EXPLICIT
         call XPLB_ABQERR(-3,'Maximum number of RMAP iterations'//
     . ' reached. Maximum number of iterations: %I, abs(f) = %R',
     .                    maxIter,abs(f),)
#else
         write(*,*) 'Maximum number of iterations: ',maxIter
         write(*,*) 'abs(f) = ',abs(f)
         error stop 'ERROR: Maximum number of RMAP iterations reached'
#endif
            endif
!-----------------------------------------------------------------------
!       Updating damage
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 1
            call UpdateDamage(VVF,sigma,dgamma,q1,q2)
            if((VVF.ge.VVFC).or.(VVF.ge.one))then
              VVF = min(VVFC,one)
              isActive = 0
            endif
#endif
          endif
!-----------------------------------------------------------------------
!       Continue if step was elastic OR end of return map
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!       Approximating the dissipated energy by using tau at iter+1
!-----------------------------------------------------------------------
          do a=1,alpha
#if SCMM_HYPO_DFLAG == 2
            Dissipation(km) = Dissipation(km)+
     .                        (one-VVF)*tau_eff(a)*dgamma(a)
#else
#if SCMM_HYPO_DFLAG == 1
            Dissipation(km) = Dissipation(km)+
     .                        (one-VVF)*tau(a)*dgamma(a)
#else
            Dissipation(km) = Dissipation(km)+tau(a)*dgamma(a)
#endif
#endif
          enddo
!-----------------------------------------------------------------------
!       Calculating incremental elastic rotation in the global coordinate system
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
! End sub-stepping
!-----------------------------------------------------------------------
        enddo! End sub-stepping
!-----------------------------------------------------------------------
!       Calculating the Cauchy stress tensor from the effective stress
!-----------------------------------------------------------------------
#if SCMM_HYPO_DFLAG == 1
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
        STATENEW(km,26)    = sqrt(half*((sigma(1)-sigma(2))**2
     +                                 +(sigma(2)-sigma(3))**2
     +                                 +(sigma(3)-sigma(1))**2)
     +                      +three*sigma(4)**2+three*sigma(5)**2
     +                      +three*sigma(6)**2)
        STATENEW(km,27) = PEQ! Equivalent von mises plastic strain
        STATENEW(km,28) = nsub! Number of sub steps
#if SCMM_HYPO_DFLAG == 1 || SCMM_HYPO_DFLAG == 2
        STATENEW(km,29) = VVF ! Damage / void volume fraction
! Is the element active or should it be deleted (Abaqus status variable)
        STATENEW(km,30) = isActive
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
      end subroutine CCCP
!-----------------------------------------------------------------------
! End preprocessor definitions
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------