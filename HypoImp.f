      include './Hypo.f'
      subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     +                RPL,DDSDDT,DRPLDE,DRPLDT,
     +                STRAN,DSTRAN,TIMEA,DTIMEA,TEMP,DTEMP,PREDEF,DPRED,
     +                CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     +                DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     +                KSPT,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
!      implicit none
      character*80 CMNAME
!      real*8 STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     +          DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     +          TIMEA(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     +          DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      real*8 SSE,SPD,SCD,PNEWDT,CELENT,RPL,DRPLDT,TEMP,DTEMP,DTIMEA
      integer NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KSTEP,
     +        KINC
!
      integer k,kk,i,j
      real*8 defgradOld(1,9),defgradNew(1,9)
      real*8 stateOld(1,NSTATV),stateNew(1,NSTATV)
      real*8 stressOld(1,6),stressNew(1,6)
      real*8 dt
      real*8 C11,C12,C44
	  real*8 Dissipation(1)! The change in dissipated inelastic specific energy (sigma_ij*D^p_ij*dt=sum(tau(alpha)*dgamma(alpha)))
	  real*8 DissipationTGT(12)
	  real*8 StressPower
	  real*8 SPRIME(6), TDROT(3,3)
!
      real*8 d(12,6),epsinc(6),spininc(3)
      real*8 L11(12),L12(12),L13(12)
      real*8 L21(12),L22(12),L23(12)
      real*8 L31(12),L32(12),L33(12)
      real*8 F(12,9),Fold(12,9)
      real*8 stateTGTold(12,NSTATV),stateTGTnew(12,NSTATV)
      real*8 stressTGTold(12,NSTATV),stressTGTnew(12,NSTATV)
      real*8 CTO(6,6),CTO2(6,6)
      real*8 pert,o2pert
      parameter(pert=1e-6,o2pert=1.0/(2.0*pert))
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
!
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
      do k=1,(NSTATV)
         stateOld(1,k) = STATEV(k)
      enddo
!-----------------------------------------------------------------------
!     Rotate the stress tensor to the right coordinate system
!-----------------------------------------------------------------------
      Call mtransp(DROT,TDROT)
	  CALL ROTSIG(STRESS,TDROT,SPRIME,1,NDI,NSHR)
	  STRESS	= SPRIME
	  stressOld(1,1) = STRESS(1)
      stressOld(1,2) = STRESS(2)
      stressOld(1,3) = STRESS(3)
      stressOld(1,4) = STRESS(4)
      stressOld(1,5) = STRESS(6)
      stressOld(1,6) = STRESS(5)
!-----------------------------------------------------------------------
!     Packaging time and time step
!-----------------------------------------------------------------------
      dt = DTIMEA
!-----------------------------------------------------------------------
!     Call the Hypo Subroutine
!-----------------------------------------------------------------------
      CALL Hypo(stressNew, stateNew, defgradNew,
     +         stressOld, stateOld, defgradOld,dt,props,
     +         1, 3, 3, NSTATV, nprops,Dissipation)
!-----------------------------------------------------------------------
!     Update Consistent tangent operator
!-----------------------------------------------------------------------
      flagCTO = int(props(17))
      if(flagCTO.eq.1)then
         C11 = props(1)
         C12 = props(2)
         C44 = props(3)
!
         DDSDDE = 0.0
!
         DDSDDE(1,1) = C11
         DDSDDE(1,2) = C12
         DDSDDE(1,3) = C12
!
         DDSDDE(2,1) = C12
         DDSDDE(2,2) = C11
         DDSDDE(2,3) = C12
!
         DDSDDE(3,1) = C12
         DDSDDE(3,2) = C12
         DDSDDE(3,3) = C11
!
         DDSDDE(4,4) = C44
         DDSDDE(5,5) = C44
         DDSDDE(6,6) = C44
!
      else
         call sinc(DFGRD0,DFGRD1,dt,epsinc,spininc)
!
         do i=1,6
            do k=1,12
               d(k,i) = epsinc(i)
            enddo
         enddo
!
         kk = 1
         do i=1,6
            do k=1,2
               d(kk,i) = epsinc(i)+pert*(-1.0)**real(k)
               kk      = kk+1
            enddo
         enddo
!
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
            F(k,1) = DFGRD0(1,1)*(L11(k)+1.0)+DFGRD0(2,1)*L12(k)
     +           +DFGRD0(3,1)*L13(k)
            F(k,4) = DFGRD0(1,2)*(L11(k)+1.0)+DFGRD0(2,2)*L12(k)
     +           +DFGRD0(3,2)*L13(k)
            F(k,9) = DFGRD0(1,3)*(L11(k)+1.0)+DFGRD0(2,3)*L12(k)
     +           +DFGRD0(3,3)*L13(k)
            F(k,7) = DFGRD0(2,1)*(L22(k)+1.0)+DFGRD0(1,1)*L21(k)
     +           +DFGRD0(3,1)*L23(k)
            F(k,2) = DFGRD0(2,2)*(L22(k)+1.0)+DFGRD0(1,2)*L21(k)
     +           +DFGRD0(3,2)*L23(k)
            F(k,5) = DFGRD0(2,3)*(L22(k)+1.0)+DFGRD0(1,3)*L21(k)
     +           +DFGRD0(3,3)*L23(k)
            F(k,6) = DFGRD0(3,1)*(L33(k)+1.0)+DFGRD0(1,1)*L31(k)
     +           +DFGRD0(2,1)*L32(k)
            F(k,8) = DFGRD0(3,2)*(L33(k)+1.0)+DFGRD0(1,2)*L31(k)
     +           +DFGRD0(2,2)*L32(k)
            F(k,3) = DFGRD0(3,3)*(L33(k)+1.0)+DFGRD0(1,3)*L31(k)
     +           +DFGRD0(2,3)*L32(k)
         enddo
!
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
!
         do i=1,NSTATV
            do k=1,12
               stateTGTold(k,i) = STATEV(i)
            enddo
         enddo
!
         do k=1,12
            stressTGTold(k,1) = STRESS(1)
            stressTGTold(k,2) = STRESS(2)
            stressTGTold(k,3) = STRESS(3)
            stressTGTold(k,4) = STRESS(4)
            stressTGTold(k,5) = STRESS(6)
            stressTGTold(k,6) = STRESS(5)
         enddo
!
         call Hypo(stressTGTnew,stateTGTnew,F,
     +         stressTGTold,stateTGTold,Fold,dt,props,
     +         12, 3, 3, NSTATV, nprops,DissipationTGT)
!
         kk = 0
         do i=1,12,2
            kk = kk+1
            do j=1,6
               CTO(kk,j) = (stressTGTnew(i+1,j)-stressTGTnew(i,j))*o2pert
            enddo
         enddo
!
         do i=1,6
            do j=1,6
               CTO2(j,i) = CTO(i,j)
            enddo
         enddo
!
         do i=1,4
            do j=1,4
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
!     Update stress and strain tensor
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
!      SPD = hisv(1,int(dm(214))+4)
!      SSE = hisv(1,int(dm(214))+1)
!-----------------------------------------------------------------------
!		Updating the dissipated inelastic specific energy
!-----------------------------------------------------------------------
	  SPD	= SPD+Dissipation(1)
!-----------------------------------------------------------------------
!		Updating the specific elastic internal energy
!-----------------------------------------------------------------------
	  StressPower = 5.d-1*(
     +          (stressOld(1,1)+stressNew(1,1))*epsinc(1) +
     +          (stressOld(1,2)+stressNew(1,2))*epsinc(2) +
     +          (stressOld(1,3)+stressNew(1,3))*epsinc(3) +
     +      2.d0*(stressOld(1,4)+stressNew(1,4))*epsinc(4) +
     +      2.d0*(stressOld(1,5)+stressNew(1,5))*epsinc(5) +
     +      2.d0*(stressOld(1,6)+stressNew(1,6))*epsinc(6) )
	  SSE	= SSE+(StressPower-Dissipation(1))
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end
