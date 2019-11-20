/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file
 * \brief Matrix type definitions
 */

#ifndef __SUPERLU_SUPERMATRIX /* allow multiple inclusions */
#define __SUPERLU_SUPERMATRIX

/********************************************
 * The matrix types are defined as follows. *
 ********************************************/
typedef enum {
    SLU_NC,    /* column-wise, no supernode */
    SLU_NCP,   /* column-wise, column-permuted, no supernode 
                  (The consecutive columns of nonzeros, after permutation,
		   may not be stored  contiguously.) */
    SLU_NR,    /* row-wize, no supernode */
    SLU_SC,    /* column-wise, supernode */
    SLU_SCP,   /* supernode, column-wise, permuted */    
    SLU_SR,    /* row-wise, supernode */
    SLU_DN,     /* Fortran style column-wise storage for dense matrix */
    SLU_NR_loc  /* distributed compressed row format  */ 
} Stype_t;

typedef enum {
    SLU_S,     /* single */
    SLU_D,     /* double */
    SLU_C,     /* single complex */
    SLU_Z      /* double complex */
} Dtype_t;

typedef enum {
    SLU_GE,    /* general */
    SLU_TRLU,  /* lower triangular, unit diagonal */
    SLU_TRUU,  /* upper triangular, unit diagonal */
    SLU_TRL,   /* lower triangular */
    SLU_TRU,   /* upper triangular */
    SLU_SYL,   /* symmetric, store lower half */
    SLU_SYU,   /* symmetric, store upper half */
    SLU_HEL,   /* Hermitian, store lower half */
    SLU_HEU    /* Hermitian, store upper half */
} Mtype_t;

typedef struct {
	Stype_t Stype; /* Storage type: interprets the storage structure 
		   	  pointed to by *Store. */
	Dtype_t Dtype; /* Data type. */
	Mtype_t Mtype; /* Matrix type: describes the mathematical property of 
			  the matrix. */
	int_t  nrow;   /* number of rows */
	int_t  ncol;   /* number of columns */
	void *Store;   /* pointer to the actual storage of the matrix */
} SuperMatrix;

/***********************************************
 * The storage schemes are defined as follows. *
 ***********************************************/

/* Stype == SLU_NC (Also known as Harwell-Boeing sparse matrix format) */
typedef struct {
    int_t  nnz;	    /* number of nonzeros in the matrix */
    void *nzval;    /* pointer to array of nonzero values, packed by column */
    int_t  *rowind; /* pointer to array of row indices of the nonzeros */
    int_t  *colptr; /* pointer to array of beginning of columns in nzval[] 
		       and rowind[]  */
                    /* Note:
		       Zero-based indexing is used;
		       colptr[] has ncol+1 entries, the last one pointing
		       beyond the last column, so that colptr[ncol] = nnz. */
} NCformat;

/* Stype == SLU_NR */
typedef struct {
    int_t  nnz;	    /* number of nonzeros in the matrix */
    void *nzval;    /* pointer to array of nonzero values, packed by raw */
    int_t  *colind; /* pointer to array of columns indices of the nonzeros */
    int_t  *rowptr; /* pointer to array of beginning of rows in nzval[] 
		       and colind[]  */
                    /* Note:
		       Zero-based indexing is used;
		       rowptr[] has nrow+1 entries, the last one pointing
		       beyond the last row, so that rowptr[nrow] = nnz. */
} NRformat;

/* Stype == SLU_SC */
typedef struct {
  int_t  nnz;	     /* number of nonzeros in the matrix */
  int_t  nsuper;     /* number of supernodes, minus 1 */
  void *nzval;       /* pointer to array of nonzero values, packed by column */
  int_t *nzval_colptr;/* pointer to array of beginning of columns in nzval[] */
  int_t *rowind;     /* pointer to array of compressed row indices of 
			rectangular supernodes */
  int_t *rowind_colptr;/* pointer to array of beginning of columns in rowind[] */
  int_t *col_to_sup;   /* col_to_sup[j] is the supernode number to which column 
			j belongs; mapping from column to supernode number. */
  int_t *sup_to_col;   /* sup_to_col[s] points to the start of the s-th 
			supernode; mapping from supernode number to column.
		        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
		              sup_to_col: 0 1 2 4 7 12           (nsuper=4) */
                     /* Note:
		        Zero-based indexing is used;
		        nzval_colptr[], rowind_colptr[], col_to_sup and
		        sup_to_col[] have ncol+1 entries, the last one
		        pointing beyond the last column.
		        For col_to_sup[], only the first ncol entries are
		        defined. For sup_to_col[], only the first nsuper+2
		        entries are defined. */
} SCformat;

/* Stype == SLU_SCP */
typedef struct {
  int_t  nnz;	     /* number of nonzeros in the matrix */
  int_t  nsuper;     /* number of supernodes */
  void *nzval;       /* pointer to array of nonzero values, packed by column */
  int_t  *nzval_colbeg;/* nzval_colbeg[j] points to beginning of column j
			  in nzval[] */
  int_t  *nzval_colend;/* nzval_colend[j] points to one past the last element
			  of column j in nzval[] */
  int_t  *rowind;      /* pointer to array of compressed row indices of 
			  rectangular supernodes */
  int_t *rowind_colbeg;/* rowind_colbeg[j] points to beginning of column j
			  in rowind[] */
  int_t *rowind_colend;/* rowind_colend[j] points to one past the last element
			  of column j in rowind[] */
  int_t *col_to_sup;   /* col_to_sup[j] is the supernode number to which column
			  j belongs; mapping from column to supernode. */
  int_t *sup_to_colbeg; /* sup_to_colbeg[s] points to the start of the s-th 
			   supernode; mapping from supernode to column.*/
  int_t *sup_to_colend; /* sup_to_colend[s] points to one past the end of the
			   s-th supernode; mapping from supernode number to
			   column.
		        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
		              sup_to_colbeg: 0 1 2 4 7              (nsuper=4)
			      sup_to_colend: 1 2 4 7 12                    */
                     /* Note:
		        Zero-based indexing is used;
		        nzval_colptr[], rowind_colptr[], col_to_sup and
		        sup_to_col[] have ncol+1 entries, the last one
		        pointing beyond the last column.         */
} SCPformat;

/* Stype == SLU_NCP */
typedef struct {
    int_t nnz;	  /* number of nonzeros in the matrix */
    void *nzval;  /* pointer to array of nonzero values, packed by column */
    int_t *rowind;/* pointer to array of row indices of the nonzeros */
		  /* Note: nzval[]/rowind[] always have the same length */
    int_t *colbeg;/* colbeg[j] points to the beginning of column j in nzval[] 
                     and rowind[]  */
    int_t *colend;/* colend[j] points to one past the last element of column
		     j in nzval[] and rowind[]  */
		  /* Note:
		     Zero-based indexing is used;
		     The consecutive columns of the nonzeros may not be 
		     contiguous in storage, because the matrix has been 
		     postmultiplied by a column permutation matrix. */
} NCPformat;

/* Stype == SLU_DN */
typedef struct {
    int_t lda;    /* leading dimension */
    void *nzval;  /* array of size lda*ncol to represent a dense matrix */
} DNformat;

/* Stype == SLU_NR_loc (Distributed Compressed Row Format) */
typedef struct {
    int_t nnz_loc;   /* number of nonzeros in the local submatrix */
    int_t m_loc;     /* number of rows local to this processor */
    int_t fst_row;   /* global index of the first row */
    void  *nzval;    /* pointer to array of nonzero values, packed by row */
    int_t *rowptr;   /* pointer to array of beginning of rows in nzval[] 
			and colind[]  */
    int_t *colind;   /* pointer to array of column indices of the nonzeros */
                     /* Note:
			Zero-based indexing is used;
			rowptr[] has n_loc + 1 entries, the last one pointing
			beyond the last row, so that rowptr[n_loc] = nnz_loc.*/
} NRformat_loc;


/* 
 *-- This contains the options used to control the solution process.
 *
 * Fact   (fact_t)
 *        Specifies whether or not the factored form of the matrix
 *        A is supplied on entry, and if not, how the matrix A should
 *        be factorizaed.
 *        = DOFACT: The matrix A will be factorized from scratch, and the
 *             factors will be stored in L and U.
 *        = SamePattern: The matrix A will be factorized assuming
 *             that a factorization of a matrix with the same sparsity
 *             pattern was performed prior to this one. Therefore, this
 *             factorization will reuse column permutation vector 
 *             ScalePermstruct->perm_c and the column elimination tree
 *             LUstruct->etree.
 *        = SamePattern_SameRowPerm: The matrix A will be factorized
 *             assuming that a factorization of a matrix with the same
 *             sparsity	pattern and similar numerical values was performed
 *             prior to this one. Therefore, this factorization will reuse
 *             both row and column scaling factors R and C, both row and
 *             column permutation vectors perm_r and perm_c, and the
 *             L & U data structures set up from the previous factorization.
 *        = FACTORED: On entry, L, U, perm_r and perm_c contain the 
 *              factored form of A. If DiagScale is not NOEQUIL, the matrix
 *              A has been equilibrated with scaling factors R and C.
 *
 * Equil  (yes_no_t)
 *        Specifies whether to equilibrate the system (scale A's row and
 *        columns to have unit norm).
 *
 * ColPerm (colperm_t)
 *        Specifies what type of column permutation to use to reduce fill.
 *        = NATURAL: use the natural ordering 
 *        = MMD_ATA: use minimum degree ordering on structure of A'*A
 *        = MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
 *        = COLAMD: use approximate minimum degree column ordering
 *        = MY_PERMC: use the ordering specified by the user
 *         
 * Trans  (trans_t)
 *        Specifies the form of the system of equations:
 *        = NOTRANS: A * X = B        (No transpose)
 *        = TRANS:   A**T * X = B     (Transpose)
 *        = CONJ:    A**H * X = B     (Transpose)
 *
 * IterRefine (IterRefine_t)
 *        Specifies whether to perform iterative refinement.
 *        = NO: no iterative refinement
 *        = SLU_SINGLE: perform iterative refinement in single precision
 *        = SLU_DOUBLE: perform iterative refinement in double precision
 *        = SLU_EXTRA: perform iterative refinement in extra precision
 *
 * DiagPivotThresh (double, in [0.0, 1.0]) (only for sequential SuperLU)
 *        Specifies the threshold used for a diagonal entry to be an
 *        acceptable pivot.
 *
 * SymmetricMode (yest_no_t)
 *        Specifies whether to use symmetric mode. Symmetric mode gives 
 *        preference to diagonal pivots, and uses an (A'+A)-based column
 *        permutation algorithm.
 *
 * PivotGrowth (yes_no_t)
 *        Specifies whether to compute the reciprocal pivot growth.
 *
 * ConditionNumber (ues_no_t)
 *        Specifies whether to compute the reciprocal condition number.
 *
 * RowPerm (rowperm_t) (only for SuperLU_DIST or ILU)
 *        Specifies whether to permute rows of the original matrix.
 *        = NO: not to permute the rows
 *        = LargeDiag: make the diagonal large relative to the off-diagonal
 *        = MY_PERMR: use the permutation given by the user
 *
 * ILU_DropRule (int)
 *        Specifies the dropping rule:
 *	  = DROP_BASIC:   Basic dropping rule, supernodal based ILUTP(tau).
 *	  = DROP_PROWS:   Supernodal based ILUTP(p,tau), p = gamma * nnz(A)/n.
 *	  = DROP_COLUMN:  Variant of ILUTP(p,tau), for j-th column,
 *			      p = gamma * nnz(A(:,j)).
 *	  = DROP_AREA:    Variation of ILUTP, for j-th column, use
 *			      nnz(F(:,1:j)) / nnz(A(:,1:j)) to control memory.
 *	  = DROP_DYNAMIC: Modify the threshold tau during factorizaion:
 *			  If nnz(L(:,1:j)) / nnz(A(:,1:j)) > gamma
 *				  tau_L(j) := MIN(tau_0, tau_L(j-1) * 2);
 *			  Otherwise
 *				  tau_L(j) := MAX(tau_0, tau_L(j-1) / 2);
 *			  tau_U(j) uses the similar rule.
 *			  NOTE: the thresholds used by L and U are separate.
 *	  = DROP_INTERP:  Compute the second dropping threshold by
 *	                  interpolation instead of sorting (default).
 *  		          In this case, the actual fill ratio is not
 *			  guaranteed to be smaller than gamma.
 *   	  Note: DROP_PROWS, DROP_COLUMN and DROP_AREA are mutually exclusive.
 *	  ( Default: DROP_BASIC | DROP_AREA )
 *
 * ILU_DropTol (double)
 *        numerical threshold for dropping.
 *
 * ILU_FillFactor (double) 
 *        Gamma in the secondary dropping.
 *
 * ILU_Norm (norm_t)
 *        Specify which norm to use to measure the row size in a
 *        supernode: infinity-norm, 1-norm, or 2-norm.
 *
 * ILU_FillTol (double)
 *        numerical threshold for zero pivot perturbation.
 *
 * ILU_MILU (milu_t)
 *        Specifies which version of MILU to use.
 *
 * ILU_MILU_Dim (double) 
 *        Dimension of the PDE if available.
 *
 * ReplaceTinyPivot (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether to replace the tiny diagonals by
 *        sqrt(epsilon)*||A|| during LU factorization.
 *
 * SolveInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        triangular solve.
 *
 * RefineInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        sparse matrix-vector multiplication routine needed in iterative
 *        refinement.
 *
 * PrintStat (yes_no_t)
 *        Specifies whether to print the solver's statistics.
 */
typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    double        DiagPivotThresh;
    yes_no_t      SymmetricMode;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    int 	  ILU_DropRule;
    double	  ILU_DropTol;    /* threshold for dropping */
    double	  ILU_FillFactor; /* gamma in the secondary dropping */
    norm_t	  ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
    double	  ILU_FillTol;    /* threshold for zero pivot perturbation */
    milu_t	  ILU_MILU;
    double	  ILU_MILU_Dim;   /* Dimension of PDE (if available) */
    yes_no_t      ParSymbFact;
    yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
    yes_no_t      PrintStat;
    int           nnzL, nnzU;      /* used to store nnzs for now       */
    int           num_lookaheads;  /* num of levels in look-ahead      */
    yes_no_t      lookahead_etree; /* use etree computed from the
				      serial symbolic factorization */
    yes_no_t      SymPattern;      /* symmetric factorization          */
} superlu_options_t;


#endif  /* __SUPERLU_SUPERMATRIX */
