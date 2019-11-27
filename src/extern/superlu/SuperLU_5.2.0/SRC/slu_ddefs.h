/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_ddefs.h
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
#ifndef __SUPERLU_dSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_dSP_DEFS

#include "superlu_config.h"
#include "slu_error.h"

/*
 * File name:		dsp_defs.h
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


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
SUPERLU_EXPORT  void
dgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
SUPERLU_EXPORT  void
dgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int *);
    /* ILU */
SUPERLU_EXPORT  void
dgsisv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
SUPERLU_EXPORT  void
dgsisx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int *);


/*! \brief Supernodal LU factor related */
SUPERLU_EXPORT  void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
dCreate_CompRow_Matrix(SuperMatrix *, int, int, int, double *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
dCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
SUPERLU_EXPORT  void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
		     Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
dCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, double *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
SUPERLU_EXPORT  void
dCopy_Dense_Matrix(int, int, double *, int, double *, int);

SUPERLU_EXPORT  void    countnz (const int, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    ilu_countnz (const int, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    fixupL (const int, const int *, GlobalLU_t *);

SUPERLU_EXPORT  void    dallocateA (int, int, double **, int **, int **);
SUPERLU_EXPORT  void    dgstrf (superlu_options_t*, SuperMatrix*,
                       int, int, int*, void *, int, int *, int *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, int *);
SUPERLU_EXPORT  int     dsnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     dsnode_bmod (const int, const int, const int, double *,
                              double *, GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  void    dpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, double *, int *, int *, int *,
			   int *, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    dpanel_bmod (const int, const int, const int, const int,
                           double *, double *, int *, int *,
			   GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  int     dcolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     dcolumn_bmod (const int, const int, double *,
			   double *, int *, int *, int,
                           GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  int     dcopy_to_ucol (int, int, int *, int *, int *,
                              double *, GlobalLU_t *);         
SUPERLU_EXPORT  int     dpivotL (const int, const double, int *, int *, 
                         int *, int *, int *, GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  void    dpruneL (const int, const int *, const int, const int,
			  const int *, const int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    dreadmt (int *, int *, int *, double **, int **, int **);
SUPERLU_EXPORT  void    dGenXtrue (int, int, double *, int);
SUPERLU_EXPORT  void    dFillRHS (trans_t, int, double *, int, SuperMatrix *,
			  SuperMatrix *);
SUPERLU_EXPORT  void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                        SuperMatrix *, SuperLUStat_t*, int *);
/* ILU */
SUPERLU_EXPORT  void    dgsitrf (superlu_options_t*, SuperMatrix*, int, int, int*,
		        void *, int, int *, int *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, int *);
SUPERLU_EXPORT  int     dldperm(int, int, int, int [], int [], double [],
                        int [],	double [], double []);
SUPERLU_EXPORT  int     ilu_dsnode_dfs (const int, const int, const int *, const int *,
			       const int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  void    ilu_dpanel_dfs (const int, const int, const int, SuperMatrix *,
			       int *, int *, double *, double *, int *, int *,
			       int *, int *, int *, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     ilu_dcolumn_dfs (const int, const int, int *, int *, int *,
				int *, int *, int *, int *, int *,
				GlobalLU_t *);
SUPERLU_EXPORT  int     ilu_dcopy_to_ucol (int, int, int *, int *, int *,
                                  double *, int, milu_t, double, int,
                                  double *, int *, GlobalLU_t *, double *);
SUPERLU_EXPORT  int     ilu_dpivotL (const int, const double, int *, int *, int, int *,
			    int *, int *, int *, double, milu_t,
                            double, GlobalLU_t *, SuperLUStat_t*);
SUPERLU_EXPORT  int     ilu_ddrop_row (superlu_options_t *, int, int, double,
                              int, int *, double *, GlobalLU_t *, 
                              double *, double *, int);


/*! \brief Driver related */

SUPERLU_EXPORT  void    dgsequ (SuperMatrix *, double *, double *, double *,
			double *, double *, int *);
SUPERLU_EXPORT  void    dlaqgs (SuperMatrix *, double *, double *, double,
                        double, double, char *);
SUPERLU_EXPORT  void    dgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, int *);
SUPERLU_EXPORT  double   dPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
SUPERLU_EXPORT  void    dgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int *, int *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, int *);

SUPERLU_EXPORT  int     sp_dtrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, double *, SuperLUStat_t*, int *);
SUPERLU_EXPORT  int     sp_dgemv (char *, double, SuperMatrix *, double *,
			int, double, double *, int);

SUPERLU_EXPORT  int     sp_dgemm (char *, char *, int, int, int, double,
			SuperMatrix *, double *, int, double, 
			double *, int);
SUPERLU_EXPORT          double dmach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
SUPERLU_EXPORT  int     dLUMemInit (fact_t, void *, int, int, int, int, int,
                            double, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, int **, double **);
SUPERLU_EXPORT  void    dSetRWork (int, int, double *, double **, double **);
SUPERLU_EXPORT  void    dLUWorkFree (int *, double *, GlobalLU_t *);
SUPERLU_EXPORT  int     dLUMemXpand (int, int, MemType, int *, GlobalLU_t *);

SUPERLU_EXPORT  double  *doubleMalloc(int);
SUPERLU_EXPORT  double  *doubleCalloc(int);
SUPERLU_EXPORT  int     dmemory_usage(const int, const int, const int, const int);
SUPERLU_EXPORT  int     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
SUPERLU_EXPORT  int     ilu_dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
SUPERLU_EXPORT  void    dreadhb(FILE *, int *, int *, int *, double **, int **, int **);
SUPERLU_EXPORT  void    dreadrb(int *, int *, int *, double **, int **, int **);
SUPERLU_EXPORT  void    dreadtriple(int *, int *, int *, double **, int **, int **);
SUPERLU_EXPORT  void    dCompRow_to_CompCol(int, int, int, double*, int*, int*,
		                   double **, int **, int **);
SUPERLU_EXPORT  void    dfill (double *, int, double);
SUPERLU_EXPORT  void    dinf_norm_error (int, SuperMatrix *, double *);
SUPERLU_EXPORT  double  dqselect(int, double *, int);


/*! \brief Routines for debugging */
SUPERLU_EXPORT  void    dPrint_CompCol_Matrix(char *, SuperMatrix *);
SUPERLU_EXPORT  void    dPrint_SuperNode_Matrix(char *, SuperMatrix *);
SUPERLU_EXPORT  void    dPrint_Dense_Matrix(char *, SuperMatrix *);
SUPERLU_EXPORT  void    dprint_lu_col(char *, int, int, int *, GlobalLU_t *);
SUPERLU_EXPORT  int     print_double_vec(char *, int, double *);
SUPERLU_EXPORT  void    dcheck_tempv(int, double *);

/*! \brief BLAS */

SUPERLU_EXPORT  int dgemm_(const char*, const char*, const int*, const int*, const int*,
                  const double*, const double*, const int*, const double*,
		  const int*, const double*, double*, const int*);
SUPERLU_EXPORT  int dtrsv_(char*, char*, char*, int*, double*, int*,
                  double*, int*);
SUPERLU_EXPORT  int dtrsm_(char*, char*, char*, char*, int*, int*,
                  double*, double*, int*, double*, int*);
SUPERLU_EXPORT  int dgemv_(char *, int *, int *, double *, double *a, int *,
                  double *, int *, double *, double *, int *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_dSP_DEFS */

