/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_sdefs.h
 * \brief Header file for real operations
 * 
 * <pre> 
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 * 
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 * </pre>
 */
#ifndef __SUPERLU_sSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_sSP_DEFS

/*
 * File name:		ssp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#endif

/* Define my integer type int_t */
typedef int int_t; /* default */

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "slu_Cnames.h"
#include "supermatrix.h"
#include "slu_util.h"

#include "superlu_config.h"

/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
SUPERLU_EXPORT  void
sgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
SUPERLU_EXPORT  void
sgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *,
       float *, float *, float *, float *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int *);
    /* ILU */
SUPERLU_EXPORT  void
sgsisv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
SUPERLU_EXPORT  void
sgsisx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *, float *, float *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int *);


/*! \brief Supernodal LU factor related */
SUPERLU_EXPORT  void
sCreate_CompCol_Matrix(SuperMatrix *, int, int, int, float *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
sCreate_CompRow_Matrix(SuperMatrix *, int, int, int, float *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
sCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
SUPERLU_EXPORT  void
sCreate_Dense_Matrix(SuperMatrix *, int, int, float *, int,
		     Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
sCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, float *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
sCopy_Dense_Matrix(int, int, float *, int, float *, int);

SUPERLU_EXPORT  void    countnz (const int, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    ilu_countnz (const int, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    fixupL (const int, const int *, GlobalLU_t *);

SUPERLU_EXPORT  void    sallocateA (int, int, float **, int **, int **);
SUPERLU_EXPORT  void    sgstrf (superlu_options_t*, SuperMatrix*,
                       int, int, int*, void *, int, int *, int *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, int *);
SUPERLU_EXPORT  int     ssnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     ssnode_bmod (const int, const int, const int, float *,
                              float *, GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  void    spanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, float *, int *, int *, int *,
			   int *, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    spanel_bmod (const int, const int, const int, const int,
                           float *, float *, int *, int *,
			   GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  int     scolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     scolumn_bmod (const int, const int, float *,
			   float *, int *, int *, int,
                           GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  int     scopy_to_ucol (int, int, int *, int *, int *,
                              float *, GlobalLU_t *);         
SUPERLU_EXPORT  int     spivotL (const int, const double, int *, int *, 
                         int *, int *, int *, GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  void    spruneL (const int, const int *, const int, const int,
			  const int *, const int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    sreadmt (int *, int *, int *, float **, int **, int **);
SUPERLU_EXPORT  void    sGenXtrue (int, int, float *, int);
SUPERLU_EXPORT  void    sFillRHS (trans_t, int, float *, int, SuperMatrix *,
			  SuperMatrix *);
SUPERLU_EXPORT  void    sgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                        SuperMatrix *, SuperLUStat_t*, int *);
/* ILU */
SUPERLU_EXPORT  void    sgsitrf (superlu_options_t*, SuperMatrix*, int, int, int*,
		        void *, int, int *, int *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, int *);
SUPERLU_EXPORT  int     sldperm(int, int, int, int [], int [], float [],
                        int [],	float [], float []);
SUPERLU_EXPORT  int     ilu_ssnode_dfs (const int, const int, const int *, const int *,
			       const int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    ilu_spanel_dfs (const int, const int, const int, SuperMatrix *,
			       int *, int *, float *, float *, int *, int *,
			       int *, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     ilu_scolumn_dfs (const int, const int, int *, int *, int *,
				int *, int *, int *, int *, int *,
				GlobalLU_t *);
SUPERLU_EXPORT  int     ilu_scopy_to_ucol (int, int, int *, int *, int *,
                                  float *, int, milu_t, double, int,
                                  float *, int *, GlobalLU_t *, float *);
SUPERLU_EXPORT  int     ilu_spivotL (const int, const double, int *, int *, int, int *,
			    int *, int *, int *, double, milu_t,
                            float, GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  int     ilu_sdrop_row (superlu_options_t *, int, int, double,
                              int, int *, double *, GlobalLU_t *, 
                              float *, float *, int);


/*! \brief Driver related */

SUPERLU_EXPORT  void    sgsequ (SuperMatrix *, float *, float *, float *,
			float *, float *, int *);
SUPERLU_EXPORT  void    slaqgs (SuperMatrix *, float *, float *, float,
                        float, float, char *);
SUPERLU_EXPORT  void    sgscon (char *, SuperMatrix *, SuperMatrix *, 
		         float, float *, SuperLUStat_t*, int *);
SUPERLU_EXPORT  float   sPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
SUPERLU_EXPORT  void    sgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int *, int *, char *, float *, 
                       float *, SuperMatrix *, SuperMatrix *,
                       float *, float *, SuperLUStat_t*, int *);

SUPERLU_EXPORT  int     sp_strsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, float *, SuperLUStat_t*, int *);
SUPERLU_EXPORT  int     sp_sgemv (char *, float, SuperMatrix *, float *,
			int, float, float *, int);

SUPERLU_EXPORT  int     sp_sgemm (char *, char *, int, int, int, float,
			SuperMatrix *, float *, int, float, 
			float *, int);
SUPERLU_EXPORT          float smach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
SUPERLU_EXPORT  int     sLUMemInit (fact_t, void *, int, int, int, int, int,
                            float, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, int **, float **);
SUPERLU_EXPORT  void    sSetRWork (int, int, float *, float **, float **);
SUPERLU_EXPORT  void    sLUWorkFree (int *, float *, GlobalLU_t *);
SUPERLU_EXPORT  int     sLUMemXpand (int, int, MemType, int *, GlobalLU_t *);

SUPERLU_EXPORT  float  *floatMalloc(int);
SUPERLU_EXPORT  float  *floatCalloc(int);
SUPERLU_EXPORT  int     smemory_usage(const int, const int, const int, const int);
SUPERLU_EXPORT  int     sQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
SUPERLU_EXPORT  int     ilu_sQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
SUPERLU_EXPORT  void    sreadhb(FILE *, int *, int *, int *, float **, int **, int **);
SUPERLU_EXPORT  void    sreadrb(int *, int *, int *, float **, int **, int **);
SUPERLU_EXPORT  void    sreadtriple(int *, int *, int *, float **, int **, int **);
SUPERLU_EXPORT  void    sCompRow_to_CompCol(int, int, int, float*, int*, int*,
		                   float **, int **, int **);
SUPERLU_EXPORT  void    sfill (float *, int, float);
SUPERLU_EXPORT  void    sinf_norm_error (int, SuperMatrix *, float *);
SUPERLU_EXPORT  float  sqselect(int, float *, int);


/*! \brief Routines for debugging */
SUPERLU_EXPORT  void    sPrint_CompCol_Matrix(char *, SuperMatrix *);
SUPERLU_EXPORT  void    sPrint_SuperNode_Matrix(char *, SuperMatrix *);
SUPERLU_EXPORT  void    sPrint_Dense_Matrix(char *, SuperMatrix *);
SUPERLU_EXPORT  void    sprint_lu_col(char *, int, int, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     print_double_vec(char *, int, double *);
SUPERLU_EXPORT  void    scheck_tempv(int, float *);

/*! \brief BLAS */

SUPERLU_EXPORT  int sgemm_(const char*, const char*, const int*, const int*, const int*,
                  const float*, const float*, const int*, const float*,
		  const int*, const float*, float*, const int*);
SUPERLU_EXPORT  int strsv_(char*, char*, char*, int*, float*, int*,
                  float*, int*);
SUPERLU_EXPORT  int strsm_(char*, char*, char*, char*, int*, int*,
                  float*, float*, int*, float*, int*);
SUPERLU_EXPORT  int sgemv_(char *, int *, int *, float *, float *a, int *,
                  float *, int *, float *, float *, int *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_sSP_DEFS */

