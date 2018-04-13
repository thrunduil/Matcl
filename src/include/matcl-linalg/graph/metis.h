/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

/// algorithm used by metis partitioner
enum class part_type
{
    edge_cut_recursive, /// minimize edge cut using multilevel recursive bisection
    edge_cut_k_way,     /// minimize edge cut using multilevel k-way partitioning
    comm_volume_k_way,  /// minimize total communication volume using multilevel
                        /// k-way partitioning
};

/// matching scheme to be used during coarsening
enum class coarsing_type
{
    random,             /// random matching
    hem                 /// sorted heavy-edge matching
};

/// algorithm used during initial partitioning
enum class initial_partitioning
{
    greedy,             /// grows a bisection using a greedy strategy
    random,             /// computes a bisection at random followed by a refinement
    edge,               /// derives a separator from an edge cut
    node                /// grow a bisection using a greedy node-based strategy
};

/// algorithm used for refinement
enum class refinement_type
{
    fm,                 /// Fiduccia-Mattheyses (FM) heuristic based cut refinement
    greedy,             /// greedy-based cut and volume refinement
    fm_2_sided,         /// two-sided node FM refinement
    fm_1_sided          /// one-sided node FM refinement
};

/// wrapper class over metis graph partitioning library (version 5.1)
class MATCL_LINALG_EXPORT metis
{
    private:
        using impl_type         = details::metis_impl;
        using impl_ptr          = std::shared_ptr<impl_type>;

        friend details::metis_impl;

    private:
        impl_ptr        m_impl;

    public:
        //------------------------------------------------------------------------
        //                      PROBLEM CONSTRUCTION
        //------------------------------------------------------------------------

        /// prepare partitioner for scalar 1.0
        metis();

        /// prepare partitioner for a symmetric matrix A; if use_weights is true,
        /// then values in the matrix will be used as weights, otherwise unit edge
        /// weights are assumed; see also set_node_weights; edge weights must be 
        /// greater than zero
        metis(const Matrix& A, bool use_weights = false);

        /// standard destructor
        ~metis();

        /// assign weights to nodes interpreted as required resources; weights are
        /// given by K x N matrix, where K is number of different resources (for example
        /// (computation power, memory) and N is number of nodes; all weights are 
        /// represented internally as integers; real weights are rescaled
        /// proportionally with scaling that makes maximum element equal to 10^3,
        /// integer weights are not changed; node weights must be greated or equal zero
        void            set_node_weights(const Matrix& W);

        /// assign communication cost associated with given node; weights are given
        /// by N x 1 matrix (N = number of nodes); all weights are represented internally
        /// as integers; real weights are rescaled proportionally with scaling that makes
        /// maximum element equal to 10^3, integer weights are not changed; node weights
        /// must be greated or equal zero
        void            set_node_size(const Matrix& W);

        /// set random number generator seed; seed is constant until is explicitly
        /// changed and therefore results are always the same
        void            set_seed(Integer s);

        /// get random number seed
        Integer         get_seed() const;

        //------------------------------------------------------------------------
        //                      COMPUTATION FUNCTIONS
        //------------------------------------------------------------------------

        /// ordering of the matrix A that reduce fill-ins for matrix factorization
        /// using multilevel nested dissection; return permutation vector p, such that
        /// A(p,p) has more sparser factorization than A
        permvec         ordering();

        /// partition a graph into k parts in order to minimize communication volume
        /// keeping workload for each part balanced; communication volume can be given
        /// by edge cut (sum of weights of edges with endpoints belonging to 
        /// different partitions) or total communication volume defined as
        ///     tot_vol = sum { s_i x nadj_i, i in Vb}
        /// where Vb is a set of nodes that have neighbores in different partition; 
        /// nadjs_i is number of partition that nodes adjacent to i belong to minus 1,
        /// ans s_i is size of i-th node (default is 1 unless is set by set_node_size);
        /// balance contraint is defined as 
        ///     max_i{ w_ji / t_ji} - 1 <= balance_j        (1)
        /// where j is constrait number, i is partition number, w_ji is a fraction of 
        /// computational resources of type j assigned to partition i defined by function
        /// set_node_weights(), t_ji is required split, and balance_j is allowed imbalance
        /// for j-th constraint;
        /// this version assumes that t_ji = 1/k (equal split) and balance_j is the same
        /// for each j and can be set by set_imbalance function. 
        Matrix          partition(Integer k, part_type obj = part_type::edge_cut_k_way);

        /// version of partition function that allows for defining target split of 
        /// computational resourses; contraints is J x k real matrix, where J is number
        /// of computational resources (determined by set_node_weights() function and 
        /// if node weights are not set, then J = 1), such that sum(set_node_weights,2) = 1;
        /// balance is J x 1 real matrix that gives imbalance tolerance for j-th resource
        /// as in (1), 0 means perfect balancing
        Matrix          partition(Integer k, const Matrix& contraints, const Matrix& balance,
                            part_type obj = part_type::edge_cut_k_way);

        /// return communitation cost obtained by last partition; -1 if partition function
        /// was not called
        Integer         objective() const;

        //------------------------------------------------------------------------
        //                          OPTIONS
        //------------------------------------------------------------------------

        /// get array that stores metis options; see metis user manual for available
        /// options; on default all options that has impact on algorithms used by meths
        /// take default values
        Integer*        get_option_ptr();
        const Integer*  get_option_ptr() const;

        /// set imbalance tolerance as in (1) used by partition function if balance
        /// tolerance is not set for each constraint; val >= 0
        metis&          set_imbalance(Real val);
        
        /// matching scheme to be used during coarsening
        metis&          set_coarsing_alg(coarsing_type type);

        /// algorithm used during initial partitioning
        metis&          set_initial_partitioning(initial_partitioning type);

        /// algorithm used for refinement.
        metis&          set_refinement_type(refinement_type type);

        /// specifies the number of different partitionings that it will compute;
        /// the final partitioning is the one that achieves the best edgecut or 
        /// communication volume; default is 1
        metis&          set_n_cuts(Integer n);

        // specifies the number of different separators that it will compute at
        /// each level of nested dissection; the final separator that is used is
        /// the smallest one; default is 1
        metis&          set_n_sep(Integer n);

        /// specifies the number of iterations for the refinement algorithms at
        /// each stage of the uncoarsening process; default is 10
        metis&          set_iter(Integer n);

        /// specifies that the partitioning routines should try to minimize the 
        /// maximum degree of the subdomain graph; if val is true then maximum
        /// connectivity is explicitly minimized
        metis&          set_min_conn(bool val);

        /// dpecifies that the coarsening will not perform any 2–hop matchings
        /// when the standard matching approach fails to sufficiently coarsen the
        /// graph; if val is true then  2–hop matching is not performed
        metis&          set_no_2hop(bool val);

        /// if val is true, then the partitioning routines should try to produce
        /// partitions that are contiguous
        metis&          set_contig(bool val);

        /// if val is true, then the graph should be compressed by combining together
        /// vertices that have identical adjacency
        metis&          set_compress(bool val);

        /// if val is true, then the connected components of the graph should first
        /// be identified and ordered separately
        metis&          set_identify_connected(bool val);

        /// nodes with degree higher than val * mean_degree will be ordered last;
        /// good values are often in the range of 6-20; val = 0 means that no nodes
        /// will be removed and ordered last
        metis&          set_pfactor(Real val);        
};

};

#pragma warning(pop)