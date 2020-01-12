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

#include "dmperm_pothen.h"
#include <memory>

int maxmatch(Integer nrows, Integer ncols, const Integer* colstr0, const Integer* rowind0, 
             Integer* prevcl, Integer* prevrw, Integer* marker, Integer* tryrow, Integer* nxtchp, 
             Integer*  rowset, Integer*  colset)
{
    ::memset(marker,0,nrows*sizeof(Integer));

    //integer        nodec, col, nextrw, lastrw, xrow, row, nxtcol,
    // $               prow, pcol

    Integer nextrw  = 0;
    Integer lastrw  = 0;
    Integer row     = 0;
    Integer nxtcol  = 0;
    Integer prow    = 0;
    Integer pcol    = 0;
    Integer col     = 0;

    for (Integer nodec = 1; nodec <= ncols; ++nodec)
    {
        //initialize node 'col' as the root of the path.

        col                 = nodec;
        prevrw[col-1]       = 0;
        prevcl[col-1]       = 0;
        nxtchp[col-1]       = colstr0[col-1] + 1;

        // main loop begins here. Each time through, try to find a
        // cheap assignment from node col.

      lab_100:    
        nextrw              = nxtchp[col-1];
        lastrw              = colstr0[col] - 1 + 1;

        if  (nextrw > 0 )
        {
            for (Integer xrow = nextrw; xrow <= lastrw; ++xrow)
            {
                row         = rowind0[xrow-1] + 1;

                if (rowset[row-1] == 0 )
                    goto lab_400;
            };

            // mark column when all adjacent rows have been
            // considered for cheap assignment.

            nxtchp[col-1]   = -1;
        };

        // Each time through, take a step forward if possible, or
        // backtrack if not .  Quit when backtracking takes us back 
        // to the beginning of the search.

        tryrow[col-1]       = colstr0[col-1] + 1;
        nextrw              = tryrow[col-1];

        if (lastrw >= nextrw)
        {
            for (Integer xrow = nextrw; xrow <= lastrw; ++xrow)
            {
                row         = rowind0[xrow-1] + 1;
                if ( marker[row-1] < nodec )
                {
                    // row is unvisited yet for this pass.
                    // take a forward step

                    tryrow[col-1]     = xrow + 1;
                    marker[row-1]     = nodec;
                    nxtcol            = rowset[row-1];

                    if  (nxtcol < 0)
                    {
                        return 1;
                    }
                    else if ( nxtcol == col )
                    {
                        return 2;
                    }
                    else if  (nxtcol > 0)
                    {
                        // the forward step led to a matched row
                        // try to extend augmenting path from
                        // the column matched by this row.

                        prevcl[nxtcol-1]    = col;
                        prevrw[nxtcol-1]    = row;
                        tryrow[nxtcol-1]    = colstr0[nxtcol-1] + 1;
                        col                 = nxtcol;

                        goto lab_100;
                    }
                    else
                    {
                        //unmatched row
                        goto lab_400;
                    };
                };
            };
        };

        // no forward step -- backtrack.
        // if we backtrack all the way, the search is done

         col        = prevcl[col-1];
         if (col > 0 )
            goto lab_100;
         else
            goto lab_600;
 
        // update the matching by alternating the matching
        // edge backward toward the root

      lab_400:

        rowset[row-1]   = col;
        prow            = prevrw[col-1];
        pcol            = prevcl[col-1];

      lab_500:
        if (pcol > 0 )
        {
            if (rowset[prow-1] != col)
                return 3;

            rowset[prow-1]  = pcol;
            col             = pcol;
            prow            = prevrw[pcol-1];
            pcol            = prevcl[pcol-1];
            goto lab_500;
        };

      lab_600:
        ;
    };

    // compute the matching from the view of column nodes

    for (row = 1; row <= nrows; ++row)
    {
	    col                 = rowset[row-1];
	    if (col > 0)
            colset[col-1]   = row;

    };

    return 0;
}

int genbtf(Integer nrows, Integer ncols, const Integer* colstr, const Integer* rowidx, const Integer* rowstr,
           const Integer* colidx, Integer* w, Integer* rnto, Integer* cnto, Integer& nhrows, Integer& nhcols, 
           Integer& hrzcmp, Integer& nsrows, Integer& sqcmpn, Integer& nvrows, Integer& nvcols, Integer& vrtcmp, 
           Integer* rcmstr, Integer* ccmstr)
{
    // initialize
    Integer vindex      = -1;
    Integer sqindx      = -2;
    Integer hindex      = -3;

    Integer cmk         = 0;
    Integer cst         = cmk + ncols;
    Integer rmk         = cst + ncols;
    Integer rst         = rmk + nrows;
    Integer rw1         = rst + nrows;
    Integer rw2         = rw1 + nrows;
    Integer rw3         = rw2 + nrows;
    Integer cw1         = rw3 + nrows;
    Integer cw2         = cw1 + ncols;
    Integer cw3         = cw2 + ncols;

    ::memset(w, 0, (cw3 + ncols)*sizeof(Integer));

    // ... algorithm consists of three stages:
    //   1.  find a maximum matching
    //   2.  find a coarse decomposition
    //   3.  find a fine decomposition

    // find the maximum matching
    int err = maxmatch(nrows, ncols , colstr, rowidx, w + cw1, w + cmk, w + rw2, w + cw2, 
                        w + cw3, w + rst, w + cst);
    if (err)
        return err;

    for (Integer i = 0; i < nrows; ++i)
        w[rmk + i]  = sqindx;

    for (Integer i = 0; i < ncols; ++i)
        w [cmk + i] = sqindx;

    // coarse partitioning -- divide the graph into three parts

    // depth-first search from unmatched columns to get the horizontal 
    // submatrix

    err = rectblk(nrows, ncols, hindex, sqindx, colstr, rowidx,  w + cst, w + rst, w + cw1, 
                w + cw2, w + cmk, w + rmk, nhrows, nhcols );
    if (err)
        return err;

    // depth-first search from unmatched rows to get the vertical submatrix

    err = rectblk(ncols, nrows, vindex, sqindx, rowstr, colidx, w + rst, w + cst, w + rw1, 
            w + rw2, w + rmk, w + cmk, nvcols, nvrows );

    if (err)
        return err;

    // the square submatrix is what is left
    Integer nscols  = ncols - nhcols - nvcols;
    nsrows          = nrows - nhrows - nvrows;

    // begin the fine partitioning and create the new to old permutation vectors

    // find connected components in the horizontal submatrix 

    if (nhcols > 0)
    {
        Integer cmbase  = 0;
        Integer rnbase  = 0;
        Integer cnbase  = 0;

        concmp(cmbase, cnbase, rnbase, hindex, ncols , nrows , nhcols, nhrows, colstr, 
               rowidx, rowstr, colidx, w + rw1, w + cw1, w + cw2, w + rw2, w + rw3, 
               w + cw3, w + rmk, w + cmk, rcmstr, ccmstr, rnto  , cnto  , hrzcmp);
    }
    else
    {
         hrzcmp = 0;
    };

    if (nsrows > 0)
    {
        // find strongly connected components in square submatrix,
        //   putting this block into block lower triangular form.

        mmc13e(nrows , ncols , nhcols, nhrows, nsrows, sqindx, hrzcmp, rowstr, colidx, 
               w + cst, w + rw1, w + rw2, w + cw1, w + cw2, w + cmk, ccmstr, rcmstr, 
               cnto,  rnto  , sqcmpn );
    }
    else
    {         
        sqcmpn  = 0;
    };

    if (nvrows > 0)
    {
        Integer cmbase  = hrzcmp + sqcmpn;
        Integer rnbase  = nhrows + nscols;
        Integer cnbase  = nhcols + nscols;

        // find connected components in vertical submatrix

        concmp( cmbase, rnbase, cnbase, vindex, nrows , ncols , nvrows, nvcols, rowstr, 
               colidx, colstr, rowidx, w + cw1, w + rw1, w + rw2, w + cw2, w + cw3, 
               w + rw3, w + cmk, w + rmk, ccmstr, rcmstr, cnto  , rnto  , vrtcmp );
    }
    else
    {     
        vrtcmp = 0;
    };

    return 0;
};

int rectblk(Integer nrows, Integer ncols, Integer marked, Integer unmrkd, const Integer* colstr0,
            const Integer* rowidx0, const Integer* colset, const Integer* rowset, Integer* prevcl, 
            Integer* tryrow, Integer* colmrk, Integer* rowmrk, Integer& nhrows, Integer& nhcols)
{
    (void)nrows;

    nhcols          = 0;
    nhrows          = 0;
    Integer col     = 0;
    Integer row     = 0;
    Integer nextcl  = 0;

    for (Integer p = 1; p <= ncols; ++p)
    {
        // find an unmatched column to start the alternating path.

        if  (colset[p-1] == 0 )
        {
            Integer fromc   = p;

            // path starts from unmatched column "fromc"
            // put fromc into horizontal set "hc"
            // indicate fromc is the root of the path.

            nhcols          = nhcols + 1;
            colmrk[fromc-1] = marked;
            tryrow[fromc-1] = colstr0[fromc-1] + 1;
            prevcl[fromc-1] = 0;
            col             =  fromc;

            // main depth-first search loop begins here.
            // Each time through take a step forward if possible
            // or backtrack if not. quit when we backtrack to the
            // beginning of the search.

            // look for a forward step from column 'col' to an unmarked row.

            lab_100:

            Integer nextrw  = tryrow[col-1];
            for (Integer xrow = nextrw; xrow <= colstr0[col + 1-1] - 1 + 1; ++xrow)
            {
                if (rowmrk[rowidx0[xrow-1]-1+1] == unmrkd)
                {
                    // take a double forward step from 'col' to 'row'
                    // and then via matching edge from 'row' to column
                    // 'nextcl'.  ('row' must be matched since 
                    // otherwise we have found an augmenting path
                    // and the maximum matching wasn't matching.)

                    tryrow[col-1]   = xrow + 1;
                    row             = rowidx0[xrow-1] + 1;
                    rowmrk[row-1]   = marked;
                    nhrows          = nhrows + 1;

                    nextcl          = rowset[row-1];

                    if  (nextcl == 0 )
                        return 4;
      
                    nhcols          = nhcols + 1;
                    colmrk[nextcl-1]= marked;
                    prevcl[nextcl-1]= col;
                    tryrow[nextcl-1]= colstr0[nextcl-1] + 1;
                    col             = nextcl;

                    goto lab_100;
                };
      
            };

            // no forward step: backtrack.  if we backtrack
            // all the way, we have completed all searchs
            // beginning at column 'p'.

            col         = prevcl[col-1];
            if  (col != 0)
               goto lab_100;

        };
    };

    return 0;
};

void concmp(Integer cmbase, Integer rnbase, Integer cnbase, Integer vindex, Integer nrows , Integer ncols , 
           Integer nvrows, Integer nvcols, const Integer* rowstr0, const Integer* colidx0, const Integer* colstr0, 
           const Integer* rowidx0, Integer* predrw, Integer* nextrw, Integer* predcl, Integer* nextcl, 
           Integer* ctab, Integer* rtab, Integer* colmrk, Integer* rowmrk, Integer* cmclad, Integer* cmrwad, 
           Integer* cnto, Integer* rnto, Integer& numcmp )
{
    // initialization
    // cn -- the number of the scanned column node
    // rn -- the number of the scanned row node

    Integer cn      = 0;
    Integer rn      = 0;
    numcmp          = 0;

    // number of vertical rows > number of vertical columns.
    // start each search for a connected component with an unmarked
    // row in the vertical block.

    Integer col     = 0;
    Integer row     = 0;

    for (Integer p = 1; p <= nrows; ++p)
    {
        if (rowmrk[p-1] == vindex)
        {
            row     = p;

            // update the value of the current working component
            // put 'row' into the new component as the root of path

            numcmp                      = numcmp + 1;
            ctab[numcmp-1]              = cnbase + 1 + cn;
            rtab[numcmp-1]              = rnbase + 1 + rn;
            cmclad[cmbase + numcmp-1]   = ctab[numcmp-1];
            cmrwad[cmbase + numcmp-1]   = rtab[numcmp-1];
            rowmrk[row-1]               = numcmp;
            rn                          = rn + 1;
            nextrw[row-1]               = rowstr0[row-1] + 1;
            predcl[row-1]               = 0;

            // from row node to col node --
            // try to find a forward step if possible
            // else backtrack

          lab_100:            

            for (Integer xcol = nextrw[row-1]; xcol <= rowstr0[row + 1-1] -1+1; ++xcol)
            {
                col = colidx0[xcol-1]+1;

                if (colmrk[col-1] == vindex)
                {
                    // forward one step :
                    // find a forward step from row 'row' to column 
                    // 'col'.  put 'col' into the current component

                    nextrw[row-1]   = xcol + 1;
                    colmrk[col-1]   = numcmp;
                    cn              = cn + 1;
                    nextcl[col-1]   = colstr0[col-1]+1;
                    predrw[col-1]   = row;

                    goto lab_300;
                };
            };

            // backward one step  (back to col node)

            nextrw[row-1]   = rowstr0[row + 1-1]+1;
            col             = predcl[row-1];

            if (col == 0)
                goto lab_500;

            // from col node to row node try to find a forward step if possible
            // else backtrack

          lab_300:
            
            for (Integer xrow = nextcl[col-1]; xrow <= colstr0[col + 1-1] - 1+1; ++ xrow)
            {
                row     = rowidx0[xrow-1]+1;

                if (rowmrk[row-1] == vindex)
                {
                    // forward one step :
                    // find a forward step from column 'col' to
                    // row 'row'.  put row into the current component

                    nextcl[col-1]   = xrow + 1;
                    rowmrk[row-1]   = numcmp;
                    rn              = rn + 1;
                    nextrw[row-1]   = rowstr0[row-1]+1;
                    predcl[row-1]   = col;

                    goto lab_100;
                };
            };

            // backward one step  (back to row node)

            nextcl[col-1]   = colstr0[col + 1-1]+1;
            row             = predrw[col-1];
            goto lab_100;
        };

      lab_500:
        ;
    };

    // generate the column and row permutations (cnto and rnto)
    // so that each component is numbered consecutively

    cmclad[cmbase + 1 + numcmp-1] = cnbase + 1 + nvcols;
    cmrwad[cmbase + 1 + numcmp-1] = rnbase + 1 + nvrows;

    for (col = 1; col <= ncols; ++col)
    {
	    Integer compn       = colmrk[col-1];

	    if (compn > 0)
        {
            cnto[ctab[compn-1]-1]   = col;
            ctab[compn-1]           = ctab[compn-1] + 1;
            colmrk[col-1]           = vindex;
        };
    };

    for (row = 1; row <= nrows; ++row)
    {
	    Integer compn   = rowmrk[row-1];

	    if ( compn > 0 )
        {
            rnto[rtab[compn-1]-1]   = row;
            rtab[compn-1]           = rtab[compn-1] + 1;
            rowmrk[row-1]           = vindex;
        };
    };

    return;
};

void mmc13e(Integer nrows , Integer ncols , Integer nhcols, Integer nhrows, Integer nscols, Integer sqindx, 
           Integer hrzcmp, const Integer*  rowstr0, const Integer*  colind0, const Integer*  colset, 
           Integer* trycol, Integer* cbegin, Integer* lowlnk, Integer* prev, Integer* colmrk, 
           Integer* ccmstr, Integer* rcmstr, Integer* cnto, Integer*  rnto, Integer& sqcmpn)
{
    //fnlpos  is the number of pairs whose positions in final ordering
    //        have been found.
    //sqcmpn  is the number of components that have been found.
    //count   is the number of pairs on the stack  (stack pointer)

    // initialization for columns in the square partition
    Integer fnlpos  = 0;
    Integer col     = 0;
    sqcmpn          = 0;
    Integer fcol    = 0;
    Integer count   = 0;
    Integer frow    = 0;
    Integer xcol    = 0;
    Integer pair    = 0;
    Integer scol    = 0;
    Integer cmpbeg  = 0;


    for (col = 1; col <= ncols; ++col)
    {
        if (colmrk[col-1] == sqindx)
            colmrk[col-1] = 0;
    };

    for (Integer j = 1; j <= nrows; ++j)
        trycol[j-1] = rowstr0[j-1]+1;

    // look for a starting pair
    for (Integer rootcl = 1; rootcl <= ncols; ++rootcl)
    {
        if  (colmrk[rootcl-1] == 0)
        {
            // put pair (rootcl, colset(rootcl)) at beginning of stack

            fcol            = rootcl;
            count           = 1;
            lowlnk[fcol-1]  = count;
            colmrk[fcol-1]  = count;
            cbegin[nscols-1]= fcol;

            // the body of this loop puts a new pair on the stack
            // or backtracks

            for (Integer passes = 1; passes <= 2*nscols - 1; ++passes)
            {
                frow        = colset[fcol-1];
               
                // have all edges leaving pair (frow,fcol) 
                // been searched?

                if (trycol[frow-1] > 0)
                {
                    // look at edges leaving from row "frow" until
                    // we find a new column "scol" that has not
                    // yet been encountered or until all edges are
                    // exhausted.                    

                    for (xcol = trycol[frow-1]; xcol <= rowstr0[frow]-1+1; ++xcol)
                    {
                        scol        = colind0[xcol-1]+1;
                        if ( colmrk[scol-1] == 0 )
                        {
                            // put new pair  (scol, colset(scol)) on the stack

                            trycol[frow-1]  = xcol + 1;
                            prev[scol-1]    = fcol;
                            fcol            = scol;
                            count           = count + 1;
                            lowlnk[fcol-1]  = count;
                            colmrk[fcol-1]  = count;
                            Integer pos     = nscols+1-count;
                            cbegin[pos-1]   = fcol;

                            goto lab_600;
                        }
                        else if (colmrk[scol-1] > 0)
                        {
                            // has scol been on stack already?  then
                            // update value of low (fcol) if necessary

                            if (lowlnk[scol-1] < lowlnk[fcol-1])
                                lowlnk[fcol-1]  = lowlnk[scol-1];
                        };
                    };

                    // there are no more edges leaving frow
                    trycol[frow-1]  = -1;
                }

                if (lowlnk[fcol-1] >= colmrk[fcol-1])
                {
                    // is  frow  the root of a block?  if so, we have 
                    // found a component.  order the nodes in this
                    // block by starting at the top of the stack and
                    // working down to the root of the block

                    sqcmpn          = sqcmpn + 1;
                    cmpbeg          = fnlpos + 1;

                    for (Integer stackp = nscols + 1 - count; stackp <= nscols; ++stackp)
                    {
                        pair            = cbegin[stackp-1];
                        fnlpos          = fnlpos + 1;
                        colmrk[pair-1]  = fnlpos;
                        count           = count-1;
                        lowlnk[pair-1]  = nscols + 1;

                        if ( pair == fcol )
                            goto lab_500;
                    };

                    // record the starting position for the new component
    
                  lab_500:

                    cbegin[sqcmpn-1]    = cmpbeg;

                    // are there any pairs left on the stack.
                    // if so, backtrack.
                    // if not, have all the pairs been ordered?

                    if  (count == 0)
                    {
                        if (fnlpos < nscols)
                            goto lab_700;
                        else
                            goto lab_800;
                    };
                };

                // backtrack to previous pair on path
                scol        = fcol;
                fcol        = prev[fcol - 1];

                if ( lowlnk[scol-1] < lowlnk[fcol-1] )
                    lowlnk[fcol-1]  = lowlnk[scol-1];

              lab_600:
                ;
            };           
        };

      lab_700:
        ;
    };

    // put permutation in the required form
  lab_800:

    for (Integer compnt = 1; compnt <= sqcmpn; ++ compnt)
    {
        ccmstr[compnt + hrzcmp-1]   = cbegin[compnt-1] + nhcols;
        rcmstr[compnt + hrzcmp-1]   = (cbegin[compnt-1] + nhcols) - (nhcols - nhrows);
    };

    ccmstr[hrzcmp + sqcmpn + 1-1]   = nhcols + nscols + 1;
    rcmstr[hrzcmp + sqcmpn + 1-1]   = nhrows + nscols + 1;

    // note that columns not in the square partition have
    // colmrk set negative.  diagonal entries in the
    // square block all correspond to matching pairs.

    for (col = 1; col <= ncols; ++col)
    {
        Integer j   = colmrk[col-1];
        if ( j > 0 )
        {
            cnto[nhcols + j-1]  = col;
            rnto[nhrows + j-1]  = colset[col-1];
            colmrk[col-1]       = sqindx;
        };
    };

    return;
};