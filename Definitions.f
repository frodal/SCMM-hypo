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
!   If SCMM_HYPO_DFLAG is 1 then the RT damage model is used,
!   If SCMM_HYPO_DFLAG is 0 then damage and fracture is turned off.
!   By default (unless SCMM_HYPO_DFLAG is 0) the RT damage model is used.
!-----------------------------------------------------------------------
!   WARNING! The subroutines and Abaqus only supports plane strain and
!   axisymmetric elements for certain crystallographic orientations.
!   This is because Abaqus only supports stress states with 
!   sigma(1,3)==sigma(3,1)==sigma(2,3)==sigma(3,2)==0 for these 
!   elements. Use at your own risk!
!-----------------------------------------------------------------------
#define SCMM_HYPO_3D_ONLY
! #define SCMM_HYPO_2D_ONLY
! #define SCMM_HYPO_VOCE_ONLY
! #define SCMM_HYPO_KALIDINDI_ONLY
! #define SCMM_HYPO_DFLAG 0
#define SCMM_HYPO_DFLAG 1
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
#ifndef SCMM_HYPO_DFLAG
#define SCMM_HYPO_DFLAG 1
#endif
!-----------------------------------------------------------------------