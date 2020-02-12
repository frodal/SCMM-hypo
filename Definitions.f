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
!   based on hflag
!-----------------------------------------------------------------------
! #define SCMM_HYPO_VOCE_ONLY
! #define SCMM_HYPO_KALIDINDI_ONLY
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   Do not edit the lines below!
!-----------------------------------------------------------------------
#ifdef SCMM_HYPO_VOCE_ONLY
#undef SCMM_HYPO_KALIDINDI_ONLY
#endif
!-----------------------------------------------------------------------