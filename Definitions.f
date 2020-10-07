!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Preprocessor definitions for the SCMM-hypo subroutines
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   Add definitions below (uncomment and comment as you see fit)
!-----------------------------------------------------------------------
!   If SCMM_HYPO_VOCE_ONLY or SCMM_HYPO_KALIDINDI_ONLY is defined,
!   the hardening model used is determined at compile time,
!   otherwise the hardening model is determined at runtime
!   based on hflag.
!   If SCMM_HYPO_3D_ONLY is defined, only solid elements are supported.
!   If SCMM_HYPO_2D_ONLY is defined, only plane strain and axisymmetric 
!   elements are supported. Otherwise the subroutines support solid,
!   plane strain and axisymmetric elements.
!-----------------------------------------------------------------------
! #define SCMM_HYPO_3D_ONLY
! #define SCMM_HYPO_2D_ONLY
! #define SCMM_HYPO_VOCE_ONLY
! #define SCMM_HYPO_KALIDINDI_ONLY
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   Do not edit the lines below!
!-----------------------------------------------------------------------
#ifdef SCMM_HYPO_3D_ONLY
#undef SCMM_HYPO_2D_ONLY
#endif
#ifdef SCMM_HYPO_VOCE_ONLY
#undef SCMM_HYPO_KALIDINDI_ONLY
#endif
!-----------------------------------------------------------------------