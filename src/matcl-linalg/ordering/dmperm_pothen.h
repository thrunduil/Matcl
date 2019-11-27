/*  ==================================================================
              maxmatch -- find maximum matching
    ==================================================================
    maxmatch uses depth-first search to find an augmenting path from
    each column node to get the maximum matching.

    Alex Pothen and Chin-Ju Fan, Penn State University, 1988
    last modifed: Alex Pothen July 1990 
    last bcs modifications:  John Lewis, Sept. 1990

    input variables :

         nrows          number of row nodes in the graph.
         ncols          number of column nodes in the graph.
         colstr, rowind adjacency structure of graph, stored by columns

    output variables :      
         rowset         describe the matching.
                            rowset (row) = col > 0 means column "col" is matched
                                                 to row "row"
                                         = -1     means "row" is an unmatched node.

        colset          describe the matching.
                            colset (col) = row > 0  means row "row" is matched to
                                                    column "col"
                                         = -1       means "col" is an unmatched node.

     Working variables :

         prevrw (ncols) pointer toward the root of the depth-first search from a 
                        column to a row.
         prevcl (ncols) pointer toward the root of the depth-first search from a 
                        column to a column. the pair (prevrw,prevcl) represent a
                        matched pair.
         marker (nrows) marker (row) <= the index of the root of the current 
                        depth-first search.  row has been visited in current pass
                        when equality holds.
         tryrow (ncols) tryrow (col) is a pointer into rowind to the next row to 
                        be explored from column col in the depth-first search.
         nxtchp (ncols) nxtchp (col) is a pointer into rowind to the next row to be
                        explored from column col for the cheap assignment.  set to 
                        -1 when all rows have been considered for cheap assignment

Licencing:
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "matcl-core/matrix/scalar_types.h"

using Integer   = matcl::Integer;

// maxmatch uses depth-first search to find an augmenting path from
// each column node to get the maximum matching.
//
//  nrows:  number of rows
//  ncols:  number of columns;
//  dc:     column indices pointer of compressed column representation
//  dr:     row indices pointer of compressed column representation
//  prevcl: pointer toward the root of the depth-first search from a
//          column to a row (ncols x 1)
//  prevrw: pointer toward the root of the depth-first search from a
//          column to a column; the pair (prevrw,prevcl) represent a
//          matched pair (ncols x 1)
//  marker: marker (row) <= the index of the root of the current depth
//          -first search; row has been visited in current pass when 
//          equality holds (nrows x 1)
//  tryrow: tryrow (col) is a pointer into rowind to the next row to be
//          explored from column col in the depth-first search (ncols x 1)
//  nxtchp: nxtchp (col) is a pointer into rowind to the next row to be
//          explored from column col for the cheap assignment; set to -1
//          when all rows have been considered for cheap assignment
//          (ncols x 1)
//  rowset: (output) row matching, length nrows
//  colset: (output) column matching, length ncols
//  return:
//      0:  success
//      1:  internal error, search reached a forbidden column
//      2:  internal error, search followed a matching edge
//      3:  internal error, pointer toward root disagrees with  matching
//      4:  internal error, max matching is wrong
int maxmatch(Integer nrows, Integer ncols, const Integer* dc, const Integer* dr, 
             Integer* prevcl, Integer* prevrw, Integer* marker, Integer* tryrow, Integer* nxtchp, 
             Integer*  rowset, Integer*  colset);

/*
* genbtf -- find the block triangular form (dulmadge-mendelson 
*           decomposition) of a general rectangular sparse matrix
* 
* created        sept. 14, 1990 (jgl)
* last modified  oct. 4, 1990 (jgl)
* 
* algorithm by alex pothen and chin-ju fan
* this code based on code from alex pothen, penn state university
* 
* input variables:
* 
*    nrows  : number of rows in matrix
*    ncols  : number of columns in matrix
*    colstr, rowidx
*           : adjacency structure of matrix, where each column of 
*             matrix is stored contiguously (column-wise representation)
*    rowstr, colidx
*           : adjacency structure of matrix, where each row of matrix is
*             stored contiguously (row-wise representation) (yes, two 
*             copies of the matrix)
* 
* temporary storage:
* 
*    w      : integer array of length 5*nrows + 5*ncols
* 
* output variables:
* 
*    rnto   : the new to old permutation array for the rows
*    cotn   : the old to new permutation array for the columns
*    nhrows, nhcols, hrzcmp 
*           : number of rows, columns and connected components in the 
*             horizontal (underdetermined) block
*    nsrows, sqcmpn 
*           : number of rows (and columns) and strong components in the square
*             (exactly determined) block
*    nvrows, nvcols, vrtcmp 
*           : number of rows, columns and connected components in the vertical
*             (overdetermined) block
*    rcmstr : index of first row in a diagonal block (component starting row)
*             where
*                  (rcmstr(1), ..., rcmstr(hrzcmp)) give the indices for the 
*                       components in the horizontal block
*                  (rcmstr(hrzcmp+1), ..., rcmstr(hrzcmp+sqcmpn)) give the indices
*                       for the components in the square block
*                  (rcmstr(hrzcmp+sqcmpn+1), ..., 
*                   rcmstr(hrzcmp+sqcmpn+vrtcmp)) give the indices for the components
*                       in the vertical block
*                   rcmstr(hrzcmp+sqcmpn+vrtcmp+1) is equal to nrows+1 for convenience
*    ccmstr : index of first column in a diagonal block (component starting column)
*              where
*                  (ccmstr(1), ..., ccmstr(hrzcmp)) give the  indices for the components
*                       in the horizontal block
*                  (ccmstr(hrzcmp+1), ..., ccmstr(hrzcmp+sqcmpn))
*                       give the indices for the components in the square block, making 
*                       this block itself in block lower triangular form
*                  (ccmstr(hrzcmp+sqcmpn+1), ..., 
*                   ccmstr(hrzcmp+sqcmpn+vrtcmp)) give the indices for the components 
*                       in the vertical block
*                   ccmstr(hrzcmp+sqcmpn+vrtcmp+1) is equal to ncols+1 for convenience
* 
*              note -- if the matrix has entirely empty rows, these rows will be placed 
*                      in the vertical block, each as a component with one row and zero
*                      columns.  similarly, entirely empty columns will appear in the 
*                      horizontal block, each as a component with no rows and one column.
* 
* efficiency note:
* ----------------
* 
*   although it is not required by this code that the number of rows be larger 
*   than the number of columns, the first phase (the matching) will be faster 
*   in this case.  thus, in cases where the number of columns is substantially
*   larger than the number of rows, it will probably be more efficient to apply
*   this algorithm to the transpose of the matrix.  since the matrix is required
*   with both row and column representations, applying the algorithm to the
*   transposed matrix can be achieved simply by interchanging appropriate parameters
*   in the call to  genbtf.
*/
int genbtf(Integer nrows, Integer ncols, const Integer* colstr, const Integer* rowidx, const Integer* rowstr,
           const Integer* colidx, Integer* w, Integer* rnto, Integer* cnto, Integer& nhrows, Integer& nhcols, 
           Integer& hrzcmp, Integer& nsrows, Integer& sqcmpn, Integer& nvrows, Integer& nvcols, Integer& vrtcmp, 
           Integer* rcmstr, Integer* ccmstr);

/*
c     ==================================================================
c     ==================================================================
c     ====  rectblk -- find rectangular portion of matrix by        ====
c     ====             depth-first search                           ====
c     ==================================================================
c     ==================================================================

c     original -- alex pothen and chin-ju fan, penn state, 1988
c     bcs modifications, john lewis, sept. 1990

c     use a depth-first serch to find all the rows and columns, which
c     can be reached via alternating paths beginning from all the
c     unmatched columns.  comments and names describe use of code
c     for finding the 'horizontal' block.  the same code is used
c     to find the vertical block by performing exactly the same
c     operations on the transpose of the matrix.
c
c     input variables:
c
c         nrows    -- number of rows
c         ncols    -- number of columns
c         marked   -- value to store in marker vectors to indicate
c                     that row/column has been reached and is
c                     therefore in the horizontal block
c         unmrkd   -- initial value of marker vectors, indicating
c                     that row or column is free to be chosen
c         colstr, 
c         rowidx   -- adjacency structure of graph
c         colset   -- maximum matching for columns
c         rowset   -- maximum matching for rows
c
c    output variables:
c
c         nhrows  -- number of rows in horizontal block
c         nhcols  -- number of columns in horizontal block 
c         rowmrk, 
c         colmrk  -- row and column marker vectors.  
c                    = unmrkd --> row/column is in neither the
c                                  horizontal or vertical block yet
c                    = marked --> row/column has been reached via
c                                 search in this routine and lies
c                                 in the horizontal block
c                    = neither --> row/column is not free for use.
c                                  it was found to lie in another
c                                  block.
c                                  
c    working variables:
c
c         tryrow -- tryrow (col) is a pointer into rowidx to the
c                   next row to be explored from col 'col' in
c                   the search.
c         prevcl -- pointer toward the root of the search from 
c                   column to column.
*/

int rectblk(Integer nrows, Integer ncols, Integer marked, Integer unmrkd, const Integer* colstr,
            const Integer* rowidx, const Integer* colset, const Integer* rowset, Integer* prevcl, 
            Integer* tryrow, Integer* colmrk, Integer* rowmrk, Integer& nhrows, Integer& nhcols);

/*
c     ==================================================================
c     ==================================================================
c     ====  concmp -- find the connected components in the          ====
c     ====            vertical (horizontal) block                   ====
c     ==================================================================
c     ==================================================================

c     original -- alex pothen and chin-ju fan, penn state, 1988
c     bcs modifications, john lewis, sept. 19, 1990

c     concmp:  find the connected components in the subgraph spanned
c              by the rows and columns in the vertical block.  the
c              same subroutine is used to find the connected
c              components in the horizontal block -- the transpose
c              of the matrix is used for that case.
c
c     input variables:
c
c         cmbase -- the number of components found in previous fine
c                   analysis of the coarse partition
c         rnbase -- the number of rows in earlier numbered partitions
c                   (0 for the horizontal block, nhrows+nsrows for
c                    the vertical partition)
c         cnbase -- the number of columns in earlier numbered partitions
c         vindex -- used to check whether the nodes belong in the
c                   vertical block
c         nrows  -- number of rows in the matrix 
c         ncols  -- number of columns in the matrix 
c         nvrows -- number of rows in the vertical block
c         nvcols -- number of columns in the vertical block
c         rowstr, colidx
c               -- the adjacency structure of the matrix using
c                  row-wise storage
c         colstr, rowidx
c               -- the adjacency structure of the matrix using
c                  column-wise storage
c
c     output variables:
c
c        numcmp  -- number of connected components
c        colmrk  -- initially,                        
c                    colmrk(i) = vindex if i belongs to vc.
c                              < 0 otherwise.
c                    during execution, 
c                    colmrk(i) = j, if i belongs to the jth component.
c                    after execution, original values restored
c        rowmrk -- initially,                        
c                    rowmrk(i) = vindex if i belongs to vr.
c                              < 0  otherwise.
c                    during execution, 
c                    rowmrk(i) = j, if i belongs to the jth component.
c                              < 0 otherwise.
c                    after execution, original values restored
c        cmclad, cmrwad 
c               -- the address (in the new ordering) of the 
c                  first column/row in each component,
c        cnto   -- the new to old mapping for the columns
c        rnto   -- the new to old mapping for the rows
c
c     working variables:
c
c        predrw, predcl
c               -- the path stack --
c                     predrw(i) = j means that we have in the path an
c                                   edge leaving from row node j to
c                                   column node i.
c                     predcl(i) = j means that we have in the path an 
c                                   edge leaving from column node j to
c                                   row node i.
c        nextcl -- nextcl(i) is index of first unsearched edge leaving
c                      from column node i.
c        nextrw -- nextrw(i) is index of first unsearched edge leaving
c                      from row node i.
c
c        ctab, rtab
c               -- temporary copy of the address (in the new ordering)
c                  of the first column/row in each component
*/
void concmp(Integer cmbase, Integer rnbase, Integer cnbase, Integer vindex, Integer nrows , Integer ncols , 
           Integer nvrows, Integer nvcols, const Integer* rowstr, const Integer* colidx, const Integer* colstr, 
           const Integer* rowidx, Integer* predrw, Integer* nextrw, Integer* predcl, Integer* nextcl, 
           Integer* ctab, Integer* rtab, Integer* colmrk, Integer* rowmrk, Integer* cmclad, Integer* cmrwad, 
           Integer* cnto, Integer* rnto, Integer& numcmp );

/*
c     ==================================================================
c     ==================================================================
c     ====  mmc13e -- lower block triangular form of square matrix  ====
c     ==================================================================
c     ==================================================================

c     mmc13e :   modified from harwell mc13e by alex pothen and
c                chin-ju fan
c     bcs modifications, john lewis, sept. 1990

c     finds the lower block triangular form of the square submatrix
c     in the general block triangular form.  the square submatrix
c     consists entirely of matched rows and columns.  therefore,
c     with each row matched to its matching column, the submatrix
c     has a nonzero diagonal, as required by duff's algorithm.
c
c     from a graph-theoretic standard, this is the same as considering
c     the subgraph induced by sr and sc, if non-matching edges
c     are directed from rows to columns, and matching edges are shrunk 
c     into single vertices, the resulting directed graph has strongly 
c     connected components.
c
c     mmc13e uses Tarjan's algorithm to find the strongly connected
c     components by depth-first search. All the pairs have been visited
c     will be labeled in the order they are visited, and associated a
c     lowlink for each vertex, stored in the stack - lowlnk.
c
c     input variables :
c
c        nrows  -- number of rows in matrix
c        ncols  -- number of columns in matrix
c        nhcols -- number of columns in horizontal (underdetermined)
c                  partition
c        nhrows -- number of rows in horizontal (underdetermined)
c                  partition
c        nscols -- number of rows and columns in square partition
c        sqindx -- index for SR and SC, for rows and columns
c                  in the square partition
c        hrzcmp -- number of components in vertical partition
c        rowstr, colind
c               -- the adjacency structure, stored by rows
c        colset -- the row matched to a column (if any)
c
c     output variables :
c
c        sqcmpn -- number of components in the square partition
c        ccmstr -- global component start vector
c        rcmstr -- global component start vector
c        cnto   -- new to old mapping for columns
c        rnto   -- new to old mapping for rows
c
c     working variables  :
c
c        trycol -- pointer to next unsearched column for this row 
c        cbegin -- is the beginning of the component.
c        colmrk -- column mark vector.
c                  on input, is negative for all columns
c                            = sqindx for columns in sc
c                  used temporarily as a stack to
c                  store the depth-first numbering for each pair.
c                  that is, is the position of pair i in the stack
c                  if it is on it, is the permuted order of pair i for
c                  those pairs whose final position has been found and
c                  is otherwise zero for columns in sc and negative
c                  for all other columns.
c                  on output, is restored to original values
c        lowlnk -- stores the lowlink for each pair.
c                  is the smallest stack position of any pair to which
c                  a path from pair i has been found. it is set to n+1
c                  when pair i is removed from the stack.
c   	 prev   -- is the pair at the end of the path when pair i was 
c                  placed on the stack
*/

void mmc13e(Integer nrows , Integer ncols , Integer nhcols, Integer nhrows, Integer nscols, Integer sqindx, 
           Integer hrzcmp, const Integer*  rowstr, const Integer*  colind, const Integer*  colset, 
           Integer* trycol, Integer* cbegin, Integer* lowlnk, Integer* prev, Integer* colmrk, 
           Integer* ccmstr, Integer* rcmstr, Integer* cnto, Integer*  rnto, Integer& sqcmpn);
