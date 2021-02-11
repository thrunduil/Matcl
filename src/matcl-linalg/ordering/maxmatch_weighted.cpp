/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "maxmatch_weighted.h"
#include <algorithm>
#include "matcl-core/lib_functions/constants.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace constants
{

template<>
Integer inf<Integer>()
{
    return matcl::constants::max_int();
};

}};

namespace matcl { namespace details
{

//minimum sum assignment problem as per Duff and Koster.
// This code is adapted from MC64 v 1.6.0 by: Jonathan Hogg (SPRAL library)
// converted to C++ by Pawel Kowal

/*
SPRAL: Sparse Parallel Robust Algorithm Library

About

The Sparse Parallel Robust Algorithm Library (SPRAL) is an open-source library for sparse 
linear algebra and related algorithms. It is primarily developed and maintained by the Numerical 
Analysis Group at STFC’s Rutherford Appleton Laboratory.

Further information can be found at the project website: http://www.numerical.rl.ac.uk/spral

Licence

At present, all code and documentation is covered by the following 3-clause BSD licence.
Future additions may use a di?erent licence.

Copyright (c) 2014, The Science and Technology Facilities Council (STFC)  
All rights reserved.  
 
Redistribution and use in source and binary forms, with or without  
modification, are permitted provided that the following conditions are met:  
    * Redistributions of source code must retain the above copyright  
      notice, this list of conditions and the following disclaimer.  
    * Redistributions in binary form must reproduce the above copyright  
      notice, this list of conditions and the following disclaimer in the  
      documentation and/or other materials provided with the distribution.  
    * Neither the name of the STFC nor the names of its contributors may be  
      used to endorse or promote products derived from this software without  
      specific prior written permission.  
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED  
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE  
DISCLAIMED. IN NO EVENT SHALL STFC BE LIABLE FOR ANY DIRECT, INDIRECT,  
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT  
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE  
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF  
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

template<class V>
void hungarian_init_heurisitic(bool min, Integer m, Integer n, const Integer* ptr, const Integer* row, 
        const V* val, Integer& num, Integer* iperm, Integer* jperm, V* dualu, V* d, Integer* work_l,
        Integer* work_search_from);

template<class V>
void heap_update(Integer idx, Integer N, Integer* Q, V* D, Integer* L);

template<class V>
Integer heap_pop(Integer& qlen, Integer m, Integer* Q, V* D, Integer* L);

template<class V>
void heap_delete(Integer lpos, Integer& qlen, Integer m, Integer* Q, V* D, Integer* L);

template<class V>
void details::hungarian_match(bool minimum, Integer m, Integer n, const Integer* ptr_c, const Integer* rows, 
                const V* vals, Integer* iperm, Integer* jperm, Integer& num, V* dualu, V* dualv, Integer* iwork,
                V* work)
{
    Integer* l      = iwork + m * 0;
    Integer* q      = iwork + m * 1;
    Integer* pr     = iwork + m * 2 + n * 0;
    Integer* out    = iwork + m * 2 + n * 1;
    Integer* search = iwork + m * 2 + n * 2;

    iwork           = iwork + m * 2 + n * 3;

    V*      d       = work;
    work            = d + std::max(m,n);

    num             = 0;
    for (Integer i = 0; i < m; ++i)
        iperm[i]    = 0;
    
    for (Integer i = 0; i < n; ++i)
        jperm[i]    = 0;

   hungarian_init_heurisitic(minimum, m, n, ptr_c, rows, vals, num, iperm, jperm, dualu, d, 
                             l, search);

   if(num == std::min(m,n)) 
       goto lab_1000; // If we got a complete matching, we're done

    // Repeatedly find augmenting paths until all columns are included in the
    // matching. At every step the current matching is optimal on the restriction
    // of the graph to currently matched rows and columns.

    // Main loop ... each pass round this loop is similar to Dijkstra's
    // algorithm for solving the single source shortest path problem
    for (Integer i = 0; i < m; ++i)
        d[i]        = constants::inf<V>();

    for (Integer i = 0; i < m; ++i)
        l[i]        = 0;
   
    Integer isp     = -1; 
    Integer jsp     = -1;

    for (Integer jord = 1; jord <= n; ++jord)
    {
        if (jperm[jord-1] != 0) 
            continue;

        // jord is next unmatched column
        // dmin is the length of shortest path in the tree
        V dmin              = constants::inf<V>();
        Integer qlen        = 0;
        Integer low         = m + 1;
        Integer up          = m + 1;

        // csp is the cost of the shortest augmenting path to unassigned row
        // row(isp). The corresponding column index is jsp.
        V csp               = constants::inf<V>();

        // Build shortest path tree starting from unassigned column (root) jord
        Integer j           = jord;
        pr[j-1]             = -1;

        // Scan column j
        for (Integer k = ptr_c[j-1]+1; k <= ptr_c[j]; ++k)
        {
            Integer i       = rows[k-1]+1;
            V vt            = minimum? vals[k-1] : -vals[k-1];
            V dnew          = vt - dualu[i-1];
            
            if(dnew >= csp) 
                continue;

            if (iperm[i-1] == 0)
            {
                csp         = dnew;
                isp         = k;
                jsp         = j;
            }
            else
            {
                if(dnew < dmin)
                    dmin    = dnew;

                d[i-1]      = dnew;
                qlen        = qlen + 1;
                q[qlen-1]   = k;
            };
        };

        // Initialize heap Q and Q2 with rows held in q(1:qlen)
        Integer q0          = qlen;
        qlen                = 0;
        
        for (Integer kk = 1; kk <= q0; ++kk)
        {
            Integer k       = q[kk-1];
            Integer i       = rows[k-1]+1;

            if (csp <= d[i-1])
            {
                d[i-1]      = constants::inf<V>();
                continue;
            };

            if (d[i-1] <= dmin)
            {
                low         = low - 1;
                q[low-1]    = i;
                l[i-1]      = low;
            }
            else
            {
                qlen        = qlen + 1;
                l[i-1]      = qlen;
                
                heap_update<V>(i,m,q,d,l);
            };

            // Update tree
            Integer jj      = iperm[i-1];
            out[jj-1]       = k;
            pr[jj-1]        = j;
        };

        for (Integer jdum = 1; jdum <= num; ++jdum)
        {
            // If Q2 is empty, extract rows from Q
            if(low == up) 
            {
                if(qlen == 0)
                    break;

                // Peek at top of heap
                Integer i   = q[1-1];

                if(d[i-1] >= csp)
                    break;

                dmin        = d[i-1];
                
                // Extract all paths that have length dmin and store in q(low:up-1)
            
                while(qlen > 0)
                {
                    // Peek at top of heap
                    i   = q[1-1];

                    if (d[i-1] > dmin)
                        break;

                    i           = heap_pop(qlen,m,q,d,l);
                    low         = low - 1;
                    q[low-1]    = i;
                    l[i-1]      = low;
                };
            };

            // q0 is row whose distance d(q0) to the root is smallest
            q0                  = q[up-1-1];
            V dq0               = d[q0-1];

            // Exit loop if path to q0 is longer than the shortest augmenting path
            if(dq0 >= csp) 
                break;

            up                  = up - 1;

            // Scan column that matches with row q0
            j                   = iperm[q0-1];
            V vt                = minimum ? vals[jperm[j-1]-1] : -vals[jperm[j-1]-1];
            V vj                = dq0 - vt + dualu[q0-1];

            for (Integer k = ptr_c[j-1]+1; k <= ptr_c[j]; ++k)
            {
                Integer i       = rows[k-1]+1;

                if(l[i-1] >= up)
                    continue;

                // dnew is new cost
                vt              = minimum ? vals[k-1] : -vals[k-1];
                V dnew          = vj + vt - dualu[i-1];

                // Do not update d[i-1] if dnew ge cost of shortest path
                if(dnew >= csp) 
                    continue;
                
                if (iperm[i-1] == 0)
                {
                    // Row i is unmatched; update shortest path info
                    csp = dnew;
                    isp = k;
                    jsp = j;
                }
                else
                {
                    // Row i is matched; do not update d[i-1] if dnew is larger
                    V di    = d[i-1];

                    if(di <= dnew)
                        continue;
                    if(l[i-1] >= low) 
                        continue;

                    d[i-1]  = dnew;
                    
                    if(dnew <= dmin)
                    {
                        Integer lpos    = l[i-1];

                        if(lpos != 0) 
                            heap_delete(lpos,qlen,m,q,d,l);

                        low             = low - 1;
                        q[low-1]        = i;
                        l[i-1]          = low;
                    }
                    else
                    {
                        if(l[i-1] == 0)
                        {
                            qlen    = qlen + 1;
                            l[i-1]  = qlen;
                        }
                        
                        // d[i-1] has changed
                        heap_update(i,m,q,d,l);
                    }

                    // Update tree
                    Integer jj  = iperm[i-1];
                    out[jj-1]   = k;
                    pr[jj-1]    = j;                        
                }                
            }
        };

        // If csp = INF, no augmenting path is found
        if(csp == constants::inf<V>()) 
            goto lab_190;

        // Find augmenting path by tracing backward in pr; update iperm,jperm
        num             = num + 1;
        Integer i       = rows[isp-1]+1;
        iperm[i-1]      = jsp;
        jperm[jsp-1]    = isp;
        j               = jsp;

        for (Integer jdum = 1; jdum <= num; ++jdum)
        {
            Integer jj  = pr[j-1];

            if(jj == -1) 
                break;

            Integer k   = out[j-1];
            i           = rows[k-1]+1;
            iperm[i-1]  = jj;
            jperm[jj-1] = k;
            j           = jj;
        };

        // Update U for rows in q(up:m)
        for (Integer kk = up; kk <= m; ++kk)
        {
            i           = q[kk-1];
            dualu[i-1]  = dualu[i-1] + d[i-1] - csp;
        };

      lab_190:
        for (Integer kk = low; kk <= m; ++kk)
        {
            i           = q[kk-1];
            d[i-1]      = constants::inf<V>();
            l[i-1]      = 0;
        };
        for (Integer kk = 1; kk <= qlen; ++kk)
        {
            i           = q[kk-1];
            d[i-1]      = constants::inf<V>();
            l[i-1]      = 0;
        };
    };

  lab_1000:

    // Set dual column variables
    for (Integer j = 1; j <= n; ++j)
    {
        Integer k       = jperm[j-1];

        V vt            = minimum ? vals[k-1] : -vals[k-1];
        if(k != 0)
            dualv[j-1]  = vt - dualu[rows[k-1]];
        else
            dualv[j-1]  = V(0.0);
    };

    // Zero dual row variables for unmatched rows
    for (Integer i = 0; i < m; ++i)
    {
        if (iperm[i] == 0)
            dualu[i]    = V(0.0);
    };

    for (Integer i = 0; i < n; ++i)
        jperm[i]    = 0;

    for (Integer i = 1; i <= m; ++i)
    {
        if (iperm[i-1] == 0)
        {
            continue;
        }
        else
        {
            Integer j   = iperm[i-1];
            jperm[j-1]  = i;
        };
    };

    if (minimum == false)
    {
        for (Integer j = 0; j < m; ++j)
            dualu[j]    = -dualu[j];
        
        for (Integer j = 0; j < n; ++j)
            dualv[j]    = -dualv[j];
    };
};

// Subroutine that initialize matching and (row) dual variable into a suitbale
// state for main Hungarian algorithm.
//
// The heuristic guaruntees that the generated partial matching is optimal
// on the restriction of the graph to the matched rows and columns.

template<class V>
void hungarian_init_heurisitic(bool minimum, Integer m, Integer n, const Integer* cols, 
        const Integer* rows, const V* vals, Integer& num, Integer* iperm, Integer* jperm, 
        V* dualu, V* work_d, Integer* work_l, Integer* work_search_from)
{
    // Set up initial matching on smallest entry in each row (as far as possible)
    //
    // Find smallest entry in each col, and record it
    for (Integer i = 0; i < m; ++i)
        dualu[i]    = constants::inf<V>();

    for (Integer i = 0; i < m; ++i)
        work_l[i]   = -1;

    V mult          = minimum ? V(1) : V(-1);

    for (Integer j = 0; j < n; ++j)
    {
        for (Integer k = cols[j]; k < cols[j + 1]; ++k)
        {
            Integer i   = rows[k];
            V vt        = mult * vals[k];
            
            if (vt > dualu[i])
                continue;

            dualu[i]    = vt;       // Initialize dual variables
            iperm[i]    = j + 1;    // Record col
            work_l[i]   = k;        // Record posn in row(:)
        };
    };

    // Loop over rows in turn. If we can match on smallest entry in row (i.e.
    // column not already matched) then do so. Avoid matching on dense columns
    // as this makes Hungarian algorithm take longer.
    for (Integer i = 0; i < m; ++i)
    {
        // Smallest entry in row i is (i,j)
        Integer j   = iperm[i] - 1;

         // skip empty rows
        if(j == -1)
            continue;

        iperm[i]    = 0;

         // If we've already matched column j, skip this row
        if (jperm[j] != 0)
            continue;

        // Don't choose cheap assignment from dense columns
        if (cols[j + 1] - cols[j]  >  m/10 && m > 50)
            continue;

        // Assignment of column j to row i
        num         = num + 1;
        iperm[i]    = j + 1;
        jperm[j]    = work_l[i] + 1;
    };

    // If we already have a complete matching, we're already done
    if (num == std::min(m,n))
        return;

    // Scan unassigned columns; improve assignment
    for (Integer i = 0; i < n; ++i)
        work_d[i]   = V(0.0);

    for (Integer i = 0; i < n; ++i)
        work_search_from[i] = cols[i];
    
    //improve_assign: &
    for (Integer j = 0; j < n; ++j)
    {
         // column j already matched
        if (jperm[j] != 0)
            continue;

        // column j is empty
        if (cols[j] >= cols[j + 1])
            continue;

        // Find smallest value of di = a_ij - u_i in column j
        // In case of a tie, prefer first unmatched, then first matched row.
        Integer i0  = rows[cols[j]];
        V vt        = mult * vals[cols[j]];
        V vj        = vt - dualu[i0];
        Integer k0  = cols[j];

        for (Integer k = cols[j] + 1; k < cols[j + 1]; ++k)
        {
            Integer i   = rows[k];
            vt          = mult * vals[k];
            V di        = vt - dualu[i];

            if(di > vj)
                continue;

            if(di == vj && di != constants::inf<V>())
            {
                if(iperm[i] != 0 || iperm[i0] == 0)
                    continue;
            };

            vj          = di;
            i0          = i;
            k0          = k;
        };

        // Record value of matching on (i0,j)
        work_d[j]       = vj;

        // If row i is unmatched, then match on (i0,j) immediately
        if (iperm[i0] == 0)
        {
            num                 = num + 1;
            jperm[j]            = k0 + 1;
            iperm[i0]           = j + 1;
            work_search_from[j] = k0 + 1;
            continue;
        };

        // Otherwise, row i is matched. Consider all rows i in column j that tie
        // for this vj value. Such a row currently matches on (i,jj). Scan column
        // jj looking for an unmatched row ii that improves value of matching. If
        // one exists, then augment along length 2 path (i,j)->(ii,jj)
        
        for (Integer k = k0; k < cols[j + 1]; ++k)
        {
            Integer i   = rows[k];
            vt          = mult * vals[k];

            // Not a tie for vj value
            if (vt - dualu[i] > vj)
                continue; 

            Integer jj = iperm[i] - 1;
            // Scan remaining part of assigned column jj
            
            for (Integer kk = work_search_from[jj]; kk < cols[jj + 1]; ++kk)
            {                
                Integer ii  = rows[kk];

                if (iperm[ii] > 0)
                    continue; // row ii already matched

                vt          = mult * vals[kk];

                if (vt - dualu[ii] <= work_d[jj])
                {
                    // By matching on (i,j) and (ii,jj) we do better than existing
                    // matching on (i,jj) alone.
                    jperm[jj]           = kk + 1;
                    iperm[ii]           = jj + 1;
                    work_search_from[jj]= kk + 1;
                    num                 = num + 1;
                    jperm[j]            = k + 1;
                    iperm[i]            = j + 1;
                    work_search_from[j] = k + 1;

                    goto lab_improve_assign;
                };
            };

            work_search_from[jj]    = cols[jj + 1];
        };

      lab_improve_assign:
        ;
    };
};

// The root node is deleted from the binary heap.
// This code is adapted from MC64 v 1.6.0

template<class V>
Integer heap_pop(Integer& QLEN, Integer N, Integer* Q, V* val, Integer* L)
{
    // Return value is the old root of the heap
    Integer ret = Q[0];

    // Delete the root
    heap_delete(1, QLEN, N, Q, val, L);
    return ret;
};

// Value associated with index i has decreased, update position in heap
// as approriate.
// This code is adapted from MC64 v 1.6.0
template<class V>
void heap_update(Integer idx, Integer N, Integer* Q, V* val, Integer* L)
{
    (void)N;

    //Get current position of i in heap
    Integer pos = L[idx-1];

    if (pos <= 1)
    {
        // idx is already at root of heap, but set q as it may have only just
        // been inserted.
        Q[pos-1]    = idx;
        return;
    };
    
    // Keep trying to move i towards root of heap until it can't go any further
    V v     = val[idx-1];

    // while not at root of heap
    while (pos > 1)
    {
        Integer parent_pos  = pos / 2;
        Integer parent_idx  = Q[parent_pos-1];

        // If parent is better than idx, stop moving
        if (v >= val[parent_idx-1]) 
            break;

        // Otherwise, swap idx and parent
        Q[pos-1]        = parent_idx;
        L[parent_idx-1] = pos;
        pos             = parent_pos;
    };

    // Finally set idx in the place it reached.
    Q[pos-1]    = idx;
    L[idx-1]    = pos;
};

// Delete element in poisition pos0 from the heap
// This code is adapted from MC64 v 1.6.0
template<class V>
void heap_delete(Integer pos0, Integer& QLEN, Integer N, Integer* Q, V* D, Integer* L)
{
    (void)N;

    // If we're trying to remove the last item, just delete it.
    if (QLEN == pos0)
    {
        QLEN    = QLEN - 1;
        return;
    };

    // Replace index in position pos0 with last item and fix heap property
    Integer idx = Q[QLEN-1];
    V v         = D[idx-1];

    // shrink heap
    QLEN        = QLEN - 1;
    // pos is current position of node I in the tree
    Integer pos = pos0;

    // Move up if appropriate
    if (pos > 1)
    {
        for(;;)
        {
            Integer parent  = pos / 2;
            Integer QK      = Q[parent-1];

            if (v >= D[QK-1])
                break;

            Q[pos-1]        = QK;
            L[QK-1]         = pos;
            pos             = parent;
            if (pos <= 1)
                break;
        };
    };

    Q[pos-1]    = idx;
    L[idx-1]    = pos;

    // Item moved up, hence doesn't need to move down
    if (pos != pos0)
        return;
    
    // Otherwise, move item down
    for (;;)
    {
        Integer child   = 2 * pos;

        if (child > QLEN)
            break;

        V DK            = D[Q[child-1]-1];

        if (child < QLEN)
        {
            V DR        = D[Q[child]-1];

            if (DK > DR)
            {
                child   = child + 1;
                DK      = DR;
            };
        };

        if (v <= DK)
            break;

        Integer QK  = Q[child-1];
        Q[pos-1]    = QK;
        L[QK-1]     = pos;
        pos         = child;
    };

    Q[pos-1]    = idx;
    L[idx-1]    = pos;
};

template
void hungarian_match(bool min, Integer m, Integer n, const Integer* ptr_c, const Integer* ptr_r, const Real* val,
                     Integer* iperm, Integer* jperm, Integer& num, Real* dualu, Real* dualv, 
                     Integer* iwork, Real* work);
template
void hungarian_match(bool min, Integer m, Integer n, const Integer* ptr_c, const Integer* ptr_r, const Float* val,
                     Integer* iperm, Integer* jperm, Integer& num, Float* dualu, Float* dualv, Integer* iwork,
                     Float* work);
template
void hungarian_match(bool min, Integer m, Integer n, const Integer* ptr_c, const Integer* ptr_r, const Integer* val,
                     Integer* iperm, Integer* jperm, Integer& num, Integer* dualu, Integer* dualv, 
                     Integer* iwork, Integer* work);

}};