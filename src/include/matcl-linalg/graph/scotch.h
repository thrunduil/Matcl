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
#pragma once

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/details/linalg_fwd.h"

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/config_linalg.h"

#pragma warning(push)
#pragma warning(disable: 4251)  //needs to have dll-interface

namespace matcl
{

/// general strategy used by scotch_partitioner
enum class scotch_strategy
{
    default_val,    /// default behavior
    balance,        /// enforce load balance as much as possible
    quality,        /// privilege quality over speed
    speed,          /// privilege speed over quality
};

/// wrapper class over scotch graph partitioning library
class MATCL_LINALG_EXPORT scotch
{
    private:
        using impl_type         = details::scotch_partit_impl;
        using impl_ptr          = std::shared_ptr<impl_type>;

        friend details::scotch_partit_impl;

    public:
        using order_ext_ret     = tuple<permvec, Integer, Matrix, Matrix>;

    private:
        impl_ptr        m_impl;

    public:
        /// prepare partitioner for scalar 1.0
        scotch();

        /// prepare partitioner for a symmetric matrix A; if use_weights is true,
        /// then values in the matrix will be used as weights, otherwise unit edge
        /// weights are assumed; see also set_node_weights
        scotch(const Matrix& A, bool use_weights = false);

        /// standard destructor
        ~scotch();

        /// set general strategy for all graph partitioning / ordering / etc. functions
        void            set_strategy(scotch_strategy strategy);

        /// assign weights to nodes; if not set then unit weights are assumed; all 
        /// weights are represented internally as integers; real weights are rescaled
        /// proportionally with scaling that makes maximum element equal to 1000,
        /// integer weights are not changed
        void            set_node_weights(const Matrix& W);

        /// set random number generator seed; seed is constant until is explicitly
        /// changed and therefore results are always the same
        void            set_seed(Integer s);

        /// get random number seed
        Integer         get_seed() const;

        /// given acyclic graph and weights assigned to nodes (interpreted as computation
        /// cost), edges (interpreted as communication cost), and processors (intepreped
        /// as computation power) assign function tries to assign nodes to K processors in
        /// order to:
        ///     1. minimalize communication cost between groups
        ///     2. keep computation costs in given group within specified constraints
        ///
        /// proc_weights: optional weights assigned to processors, if not given then unit 
        /// weights are assumed;
        /// load_imbalance: optional allowed mean absolute difference between assumed workload 
        /// share on processors and obtained by partitioner (value in [0,1])
        /// return N x K aggregation matrix C; some groups may contain zero elements
        Matrix          assign(Integer K, Real load_balande = 0.01);
        Matrix          assign(Integer K, const Matrix& proc_weights, Real load_imbalance = 0.01);

        /// assign with preassigned noded; fixed is N x 1 vector where positive value k 
        /// means, that given node is preassigned to group k, and zero or negative values
        /// means, that given node will be assigned to processor by the function
        Matrix          assign_fixed(Integer K, const Matrix& fixed, Real load_imbalance = 0.01);
        Matrix          assign_fixed(Integer K, const Matrix& proc_weights, const Matrix& fixed,
                            Real load_imbalance = 0.01);

        /// computes a coloring of the graph vertices using a variant of a Luby’s algorithm;
        /// return aggregation matrix with number of columns equal to number of colors
        Matrix          coloring();

        /// matrix is partitioned in K parts using node separators; obtained partition
        /// reoders the matrix into bordered block diagonal form; return N x (1+L) 
        /// aggregation matrix, where the first column (possibly empty) contains separator
        /// vertices and next columns describe partition into L groups (0 <= L <= K) and
        /// contain at least one element
        /// imbalance: optional parameter controlling balancing of sizes of diagonal blocks
        /// (value in [0,1])
        Matrix          node_separator(Integer K, Real imbalance = 0.05);

        /// matrix is partitioned in K parts using edge separators; return N x L 
        /// aggregation matrix, 0 <= L <= K, each group contains at least one element
        /// imbalance: optional parameter controlling balancing of sizes of blocks 
        /// (value in [0,1])
        Matrix          edge_separator(Integer K, Real imbalance = 0.05);

        /// matrix is partitioned in K parts using edge separators and some nodes can be
        /// preassigned; fixed is N x 1 vector where positive value k means, that given
        /// node is preassigned to group k, and zero or negative values means, that given
        /// node will be assigned to processor by the function; return N x L aggregation
        /// matrix, 0 <= L <= K, each group contains at least one element
        /// imbalance: optional parameter controlling balancing of sizes of blocks
        /// (value in [0,1])
        Matrix          edge_separator_fixed(Integer K, const Matrix& fixed, Real imbalance = 0.05);

        /// try to cluster nodes into k sets; return N x L aggregation matrix, 
        /// 0 <= L <= K, each group contains at least one element;
        /// imbalance: optional parameter controlling balancing of sizes of blocks
        /// (value in [0,1])
        Matrix          clustering(Integer k, Real imbalance = 0.2);

        /// ordering of the matrix A that reduce fill-ins for Cholesky factorization; 
        /// return permutation vector p, such that A(p,p) has more sparser factorization
        /// than A
        permvec         ordering();

        /// return ordering computed by ordering() function and additional information;
        /// [p, nblocks, blocks, tree] = ordering_ext()
        /// p       : permutation vector
        /// nblocks : number of column blocks in the block ordering
        /// blocks  : matrix of size nblocks + 1 describing blocks of nodes, i-th block
        ///           contain nodes with indices [blocks(i), blocks(i+1) ) (1-based)
        /// tree    : matrix of size nblocks describing separators tree; i-th element
        ///           indicates, that i-th block of nodes is subdivision of tree(i)-th
        ///           block or is a root of division if tree(i) = -1 (1-based)
        order_ext_ret   ordering_ext();
};

};

#pragma warning(pop)