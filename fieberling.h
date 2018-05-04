/*
** svn $Id: upwelling.h 795 2016-05-11 01:42:43Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Upwelling Test.
**
** Application flag:   FIBERLING
** Input script:       ocean_fiberling.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define DJ_GRADPS
#define TS_DIF2
#define MIX_S_TS
#define NO_LBC_ATT

#undef SALINITY
#define SOLVE3D
!#define AVERAGES
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV

#undef ANA_GRID
#undef ANA_INITIAL
#undef ANA_FSOBC
#undef ANA_M2OBC
#undef ANA_M3OBC
#undef ANA_TOBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SRFLUX
#define ANA_SPFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX

#undef LMD_MIXING

#ifdef LMD_MIXING
#define LMD_BKPP
#define LMD_SKPP
#define LMD_CONVEC
#define LMD_RIMIX
#define RI_SPLINES
#define LMD_SHAPIRO
#endif

#define GLS_MIXING

#ifdef GLS_MIXING
#define N2S2_HORAVG
#define K_C4ADVECTION
#define KANTHA_CLAYSON
#endif

#define SPONGE

#undef ISLAND
#ifdef ISLAND
#define ANA_INITIAL
#define ANA_GRID
#define ANA_MASK
#define MASKING
#undef QDRAG
#define LDRAG
#endif

