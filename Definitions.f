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
!   If SCMM_HYPO_DFLAG is 2 then the Han/Khadyko et al model is used,
!   If SCMM_HYPO_DFLAG is 1 then the RT damage model is used,
!   If SCMM_HYPO_DFLAG is 0 then damage and fracture is turned off.
!   By default (unless SCMM_HYPO_DFLAG is 0 or 2) the RT damage model is used.
!   If SCMM_HYPO_MODEL is 1 then the rate-dependent model is used,
!   If SCMM_HYPO_MODEL is 2 then the CCCP model is used,
!   If SCMM_HYPO_MODEL is 3 then the FC-Taylor homogenization approach
!   is used with the rate-dependent model, and
!   If SCMM_HYPO_MODEL is 4 then the FC-Taylor homogenization approach
!   is used with the CCCP model.
!   By default the rate-dependent model is used.
!   If SCMM_HYPO_NLFLAG is defined, the non-local damage model is used.
!   The non-local damage model is only supported with the rate-dependent 
!   CP model and the RT damage model in Abaqus/Explicit.
!   By default the non-local damage model is turned off.
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
! #define SCMM_HYPO_DFLAG 2
#define SCMM_HYPO_MODEL 1
! #define SCMM_HYPO_MODEL 2
! #define SCMM_HYPO_MODEL 3
! #define SCMM_HYPO_MODEL 4
#define SCMM_HYPO_NLFLAG
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
#ifndef SCMM_HYPO_MODEL
#define SCMM_HYPO_MODEL 1
#endif
#if (SCMM_HYPO_MODEL == 1 || SCMM_HYPO_MODEL == 3)  && SCMM_HYPO_DFLAG == 2
#define SCMM_HYPO_DFLAG 1
#warning "SCMM_HYPO_DFLAG == 2 is not supported with this SCMM_HYPO_MODEL"
#warning "Using SCMM_HYPO_DFLAG == 1 instead"
#endif
#if SCMM_HYPO_DFLAG != 0 && SCMM_HYPO_DFLAG != 1  && SCMM_HYPO_DFLAG != 2
#define SCMM_HYPO_DFLAG 1
#warning "The SCMM_HYPO_DFLAG used is not supported"
#warning "Using SCMM_HYPO_DFLAG == 1 instead"
#endif
#if SCMM_HYPO_MODEL != 1 && SCMM_HYPO_MODEL != 2 && SCMM_HYPO_MODEL != 3 && SCMM_HYPO_MODEL != 4
#define SCMM_HYPO_MODEL 1
#warning "The SCMM_HYPO_MODEL used is not supported"
#warning "Using SCMM_HYPO_MODEL == 1 instead"
#endif
#if SCMM_HYPO_MODEL == 1 && SCMM_HYPO_DFLAG == 2
#define SCMM_HYPO_DFLAG 1
#warning "SCMM_HYPO_DFLAG == 2 is not supported with SCMM_HYPO_MODEL == 1"
#warning "Using SCMM_HYPO_DFLAG == 1 instead"
#endif
#if defined SCMM_HYPO_NLFLAG && SCMM_HYPO_MODEL !=1 && SCMM_HYPO_DFLAG != 1
#undef SCMM_HYPO_NLFLAG
#warning "SCMM_HYPO_NLFLAG is only supported with SCMM_HYPO_MODEL == 1 and SCMM_HYPO_DFLAG == 1"
#warning "SCMM_HYPO_NLFLAG is turned off"
#endif
#if SCMM_HYPO_DFLAG != 0
#ifdef SCMM_HYPO_NLFLAG
#define SCMM_HYPO_NSTATEV 32
#else
#define SCMM_HYPO_NSTATEV 30
#endif
#else
#define SCMM_HYPO_NSTATEV 28
#endif
!-----------------------------------------------------------------------