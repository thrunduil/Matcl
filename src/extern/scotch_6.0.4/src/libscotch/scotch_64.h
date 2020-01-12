/* Copyright 2004,2007-2012,2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : library.h                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Jun-Ho HER (v6.0)                       **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : Declaration file for the LibScotch      **/
/**                static mapping and sparse matrix block  **/
/**                ordering library.                       **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     22 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to     31 may 1999     **/
/**                # Version 3.4  : from : 10 oct 1999     **/
/**                                 to     15 nov 2001     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to     20 dec 2005     **/
/**                # Version 5.0  : from : 26 apr 2006     **/
/**                                 to   : 20 feb 2008     **/
/**                # Version 5.1  : from : 30 nov 2007     **/
/**                                 to   : 07 aug 2011     **/
/**                # Version 6.0  : from : 12 sep 2008     **/
/**                                 to     28 feb 2015     **/
/**                                                        **/
/************************************************************/

#ifndef SCOTCH64_H
#define SCOTCH64_H

#include "scotch_config.h"

/*
**  The type and structure definitions.
*/

/*+ Version flags. +*/

#define SCOTCH_VERSION 6
#define SCOTCH_RELEASE 0
#define SCOTCH_PATCHLEVEL 4

/*+ Integer type. +*/

typedef int SCOTCH_Idx;

typedef int SCOTCH_Num;

#define SCOTCH_NUMMAX               ((int) (((unsigned int) 1 << ((sizeof (int) << 3) - 1)) - 1))
#define SCOTCH_NUMSTRING            "%d"

/*+ Coarsening flags +*/

#define SCOTCH_COARSENNONE          0x0000
#define SCOTCH_COARSENFOLD          0x0100
#define SCOTCH_COARSENFOLDDUP       0x0300
#define SCOTCH_COARSENNOMERGE       0x4000

/*+ Strategy string parametrization values +*/

#define SCOTCH_STRATDEFAULT         0x0000
#define SCOTCH_STRATQUALITY         0x0001
#define SCOTCH_STRATSPEED           0x0002
#define SCOTCH_STRATBALANCE         0x0004
#define SCOTCH_STRATSAFETY          0x0008
#define SCOTCH_STRATSCALABILITY     0x0010
#define SCOTCH_STRATRECURSIVE       0x0100
#define SCOTCH_STRATREMAP           0x0200
#define SCOTCH_STRATLEVELMAX        0x1000
#define SCOTCH_STRATLEVELMIN        0x2000
#define SCOTCH_STRATLEAFSIMPLE      0x4000
#define SCOTCH_STRATSEPASIMPLE      0x8000

/*+ Opaque objects. The dummy sizes of these
objects, computed at compile-time by program
"dummysizes", are given as double values for
proper padding                               +*/

typedef struct {
  double                    dummy[8];
} SCOTCH_Arch;

typedef struct {
  double                    dummy[2];
} SCOTCH_Geom;

typedef struct {
  double                    dummy[13];
} SCOTCH_Graph;

typedef struct {
  double                    dummy[15];
} SCOTCH_Mesh;

typedef struct {
  double                    dummy[4];
} SCOTCH_Mapping;

typedef struct {
  double                    dummy[12];
} SCOTCH_Ordering;

typedef struct {
  double                    dummy[1];
} SCOTCH_Strat;

/*
**  The function prototypes.
*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

SCOTCH_EXPORT SCOTCH_Arch *               SCOTCH_archAlloc    (void);
SCOTCH_EXPORT int                         SCOTCH_archInit     (SCOTCH_Arch * const);
SCOTCH_EXPORT void                        SCOTCH_archExit     (SCOTCH_Arch * const);
SCOTCH_EXPORT int                         SCOTCH_archLoad     (SCOTCH_Arch * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_archSave     (const SCOTCH_Arch * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_archBuild    (SCOTCH_Arch * const, const SCOTCH_Graph * const, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Strat * const);
SCOTCH_EXPORT char *                      SCOTCH_archName     (const SCOTCH_Arch * const);
SCOTCH_EXPORT SCOTCH_Num                  SCOTCH_archSize     (const SCOTCH_Arch * const);
SCOTCH_EXPORT int                         SCOTCH_archVar      (const SCOTCH_Arch * const);
SCOTCH_EXPORT int                         SCOTCH_archCmplt    (SCOTCH_Arch * const, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_archCmpltw   (SCOTCH_Arch * const, const SCOTCH_Num, const SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_archHcub     (SCOTCH_Arch * const, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_archMesh2    (SCOTCH_Arch * const, const SCOTCH_Num, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_archMesh3    (SCOTCH_Arch * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_archTleaf    (SCOTCH_Arch * const, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_archTorus2   (SCOTCH_Arch * const, const SCOTCH_Num, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_archTorus3   (SCOTCH_Arch * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_archVcmplt   (SCOTCH_Arch * const);
SCOTCH_EXPORT int                         SCOTCH_archVhcub    (SCOTCH_Arch * const);

SCOTCH_EXPORT void                        SCOTCH_errorProg    (const char * const);
SCOTCH_EXPORT void                        SCOTCH_errorPrint   (const char * const, ...);
SCOTCH_EXPORT void                        SCOTCH_errorPrintW  (const char * const, ...);

SCOTCH_EXPORT SCOTCH_Geom *               SCOTCH_geomAlloc    (void);
SCOTCH_EXPORT int                         SCOTCH_geomInit     (SCOTCH_Geom * const);
SCOTCH_EXPORT void                        SCOTCH_geomExit     (SCOTCH_Geom * const);
SCOTCH_EXPORT void                        SCOTCH_geomData     (const SCOTCH_Geom * const, SCOTCH_Num * const, double ** const);

SCOTCH_EXPORT SCOTCH_Graph *              SCOTCH_graphAlloc   (void);
SCOTCH_EXPORT int                         SCOTCH_graphInit    (SCOTCH_Graph * const);
SCOTCH_EXPORT void                        SCOTCH_graphExit    (SCOTCH_Graph * const);
SCOTCH_EXPORT void                        SCOTCH_graphFree    (SCOTCH_Graph * const);
SCOTCH_EXPORT int                         SCOTCH_graphLoad    (SCOTCH_Graph * const, FILE * const, const SCOTCH_Num, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_graphSave    (const SCOTCH_Graph * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphBuild   (SCOTCH_Graph * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphCoarsen (const SCOTCH_Graph * const, SCOTCH_Graph * const, SCOTCH_Num * const, const SCOTCH_Num, const double);
SCOTCH_EXPORT int                         SCOTCH_graphCoarsenBuild (const SCOTCH_Graph * const, SCOTCH_Graph * const, SCOTCH_Num * const, const SCOTCH_Num, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphColor   (const SCOTCH_Graph * const, SCOTCH_Num * const, SCOTCH_Num * const, const SCOTCH_Num);
SCOTCH_EXPORT SCOTCH_Num                  SCOTCH_graphBase    (SCOTCH_Graph * const, const SCOTCH_Num baseval);
SCOTCH_EXPORT int                         SCOTCH_graphCheck   (const SCOTCH_Graph * const);
SCOTCH_EXPORT void                        SCOTCH_graphSize    (const SCOTCH_Graph * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT void                        SCOTCH_graphData    (const SCOTCH_Graph * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num * const, SCOTCH_Num ** const, SCOTCH_Num ** const);
SCOTCH_EXPORT void                        SCOTCH_graphStat    (const SCOTCH_Graph * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const);
SCOTCH_EXPORT int                         SCOTCH_graphGeomLoadChac (SCOTCH_Graph * const, SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_graphGeomLoadHabo (SCOTCH_Graph * const, SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_graphGeomLoadMmkt (SCOTCH_Graph * const, SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_graphGeomLoadScot (SCOTCH_Graph * const, SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_graphGeomSaveChac (const SCOTCH_Graph * const, const SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_graphGeomSaveMmkt (const SCOTCH_Graph * const, const SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_graphGeomSaveScot (const SCOTCH_Graph * const, const SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);

SCOTCH_EXPORT int                         SCOTCH_graphMapInit (const SCOTCH_Graph * const, SCOTCH_Mapping * const, const SCOTCH_Arch * const, SCOTCH_Num * const);
SCOTCH_EXPORT void                        SCOTCH_graphMapExit (const SCOTCH_Graph * const, SCOTCH_Mapping * const);
SCOTCH_EXPORT int                         SCOTCH_graphMapLoad (const SCOTCH_Graph * const, const SCOTCH_Mapping * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphTabLoad (const SCOTCH_Graph * const, SCOTCH_Num * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphMapSave (const SCOTCH_Graph * const, const SCOTCH_Mapping * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphMapView (const SCOTCH_Graph * const, const SCOTCH_Mapping * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphRemapView (const SCOTCH_Graph * const, const SCOTCH_Mapping * const, const SCOTCH_Mapping * const, const double, SCOTCH_Num *, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphRemapViewRaw (const SCOTCH_Graph * const, const SCOTCH_Mapping * const, const SCOTCH_Mapping * const, const double, SCOTCH_Num *, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphMapCompute (SCOTCH_Graph * const, SCOTCH_Mapping * const, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_graphMapFixedCompute (SCOTCH_Graph * const, SCOTCH_Mapping * const, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_graphRemapCompute (SCOTCH_Graph * const, SCOTCH_Mapping * const, SCOTCH_Mapping * const, const double, const SCOTCH_Num *, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_graphRemapFixedCompute (SCOTCH_Graph * const, SCOTCH_Mapping * const, SCOTCH_Mapping * const, const double, const SCOTCH_Num *, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_graphMap     (SCOTCH_Graph * const, const SCOTCH_Arch * const, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphMapFixed (SCOTCH_Graph * const, const SCOTCH_Arch * const, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphRemap   (SCOTCH_Graph * const, const SCOTCH_Arch * const, SCOTCH_Num *, const double, const SCOTCH_Num *, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphRemapFixed (SCOTCH_Graph * const, const SCOTCH_Arch * const, SCOTCH_Num *, const double, const SCOTCH_Num *, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphPart    (SCOTCH_Graph * const, const SCOTCH_Num, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphPartFixed (SCOTCH_Graph * const, const SCOTCH_Num, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphPartOvl (SCOTCH_Graph * const, const SCOTCH_Num, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphRepart  (SCOTCH_Graph * const, const SCOTCH_Num, SCOTCH_Num * const, const double, const SCOTCH_Num *, SCOTCH_Strat * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphRepartFixed (SCOTCH_Graph * const, const SCOTCH_Num, SCOTCH_Num * const, const double, const SCOTCH_Num *, SCOTCH_Strat * const, SCOTCH_Num * const);

SCOTCH_EXPORT int                         SCOTCH_graphOrderInit (const SCOTCH_Graph * const, SCOTCH_Ordering * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT void                        SCOTCH_graphOrderExit (const SCOTCH_Graph * const, SCOTCH_Ordering * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderLoad (const SCOTCH_Graph * const, SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderSave (const SCOTCH_Graph * const, const SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderSaveMap (const SCOTCH_Graph * const, const SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderSaveTree (const SCOTCH_Graph * const, const SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderCompute (SCOTCH_Graph * const, SCOTCH_Ordering * const, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderComputeList (SCOTCH_Graph * const, SCOTCH_Ordering * const, const SCOTCH_Num, const SCOTCH_Num * const, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderFactor (const SCOTCH_Graph * const, const SCOTCH_Ordering * const, SCOTCH_Graph * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderView (const SCOTCH_Graph * const, const SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrder   (SCOTCH_Graph * const, SCOTCH_Strat * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderList (SCOTCH_Graph * const, const SCOTCH_Num, const SCOTCH_Num * const, SCOTCH_Strat * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_graphOrderCheck (const SCOTCH_Graph * const, const SCOTCH_Ordering * const);

SCOTCH_EXPORT SCOTCH_Mapping *            SCOTCH_mapAlloc     (void);

SCOTCH_EXPORT void                        SCOTCH_memFree      (void * const);
SCOTCH_EXPORT SCOTCH_Idx                  SCOTCH_memCur       (void);
SCOTCH_EXPORT SCOTCH_Idx                  SCOTCH_memMax       (void);

SCOTCH_EXPORT int                         SCOTCH_meshInit     (SCOTCH_Mesh * const);
SCOTCH_EXPORT void                        SCOTCH_meshExit     (SCOTCH_Mesh * const);
SCOTCH_EXPORT int                         SCOTCH_meshLoad     (SCOTCH_Mesh * const, FILE * const, const SCOTCH_Num);
SCOTCH_EXPORT int                         SCOTCH_meshSave     (const SCOTCH_Mesh * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_meshBuild    (SCOTCH_Mesh * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num, const SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_meshCheck    (const SCOTCH_Mesh * const);
SCOTCH_EXPORT void                        SCOTCH_meshSize     (const SCOTCH_Mesh * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT void                        SCOTCH_meshData     (const SCOTCH_Mesh * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num * const, SCOTCH_Num ** const, SCOTCH_Num * const);
SCOTCH_EXPORT void                        SCOTCH_meshStat     (const SCOTCH_Mesh * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const);
SCOTCH_EXPORT int                         SCOTCH_meshGraph    (const SCOTCH_Mesh * const, SCOTCH_Graph * const);
SCOTCH_EXPORT int                         SCOTCH_meshGeomLoadHabo (SCOTCH_Mesh * const, SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_meshGeomLoadScot (SCOTCH_Mesh * const, SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_meshGeomSaveScot (const SCOTCH_Mesh * const, const SCOTCH_Geom * const, FILE * const, FILE * const, const char * const);

SCOTCH_EXPORT int                         SCOTCH_meshOrderInit (const SCOTCH_Mesh * const, SCOTCH_Ordering * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT void                        SCOTCH_meshOrderExit (const SCOTCH_Mesh * const, SCOTCH_Ordering * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrderSave (const SCOTCH_Mesh * const, const SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrderSaveMap (const SCOTCH_Mesh * const, const SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrderSaveTree (const SCOTCH_Mesh * const, const SCOTCH_Ordering * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrderCompute (SCOTCH_Mesh * const, SCOTCH_Ordering * const, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrderComputeList (SCOTCH_Mesh * const, SCOTCH_Ordering * const, const SCOTCH_Num, const SCOTCH_Num * const, SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrder    (SCOTCH_Mesh * const, SCOTCH_Strat * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrderList (SCOTCH_Mesh * const, const SCOTCH_Num, const SCOTCH_Num * const, SCOTCH_Strat * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
SCOTCH_EXPORT int                         SCOTCH_meshOrderCheck (const SCOTCH_Mesh * const, const SCOTCH_Ordering * const);

SCOTCH_EXPORT int                         SCOTCH_numSizeof    (void);

SCOTCH_EXPORT SCOTCH_Ordering *           SCOTCH_orderAlloc   (void);

SCOTCH_EXPORT void                        SCOTCH_randomReset  (void);
SCOTCH_EXPORT void                        SCOTCH_randomSeed   (SCOTCH_Num);

SCOTCH_EXPORT SCOTCH_Strat *              SCOTCH_stratAlloc   (void);
SCOTCH_EXPORT int                         SCOTCH_stratInit    (SCOTCH_Strat * const);
SCOTCH_EXPORT void                        SCOTCH_stratExit    (SCOTCH_Strat * const);
SCOTCH_EXPORT void                        SCOTCH_stratFree    (SCOTCH_Strat * const);
SCOTCH_EXPORT int                         SCOTCH_stratSave    (const SCOTCH_Strat * const, FILE * const);
SCOTCH_EXPORT int                         SCOTCH_stratGraphBipart (SCOTCH_Strat * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_stratGraphMap (SCOTCH_Strat * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_stratGraphMapBuild (SCOTCH_Strat * const, const SCOTCH_Num, const SCOTCH_Num, const double);
SCOTCH_EXPORT int                         SCOTCH_stratGraphClusterBuild (SCOTCH_Strat * const, const SCOTCH_Num, const SCOTCH_Num, const double, const double);
SCOTCH_EXPORT int                         SCOTCH_stratGraphPartOvl (SCOTCH_Strat * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_stratGraphPartOvlBuild (SCOTCH_Strat * const, const SCOTCH_Num, const SCOTCH_Num, const double);
SCOTCH_EXPORT int                         SCOTCH_stratGraphOrder (SCOTCH_Strat * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_stratGraphOrderBuild (SCOTCH_Strat * const, const SCOTCH_Num, const SCOTCH_Num, const double);
SCOTCH_EXPORT int                         SCOTCH_stratMeshOrder (SCOTCH_Strat * const, const char * const);
SCOTCH_EXPORT int                         SCOTCH_stratMeshOrderBuild (SCOTCH_Strat * const, const SCOTCH_Num, const double);

SCOTCH_EXPORT void                        SCOTCH_version      (int * const, int * const, int * const);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* SCOTCH64_H */
