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

#include "matcl-linalg/graph/scotch.h"
//#include "matcl-linalg/petsc_utils/petsc_algs.h"
#include "matcl-linalg/general/linalg_exception.h"

//#include "petsc/mpiuni/mpi.h"

extern "C"
{
    //ptscotch.h
    #include "extern/scotch_6.0.4/include/scotch.h"
};

namespace matcl { namespace details
{

class scotch_wrapper
{
    private:
        using Mat_I     = raw::Matrix<Integer, struct_dense>;
        using arch_ptr  = std::shared_ptr<SCOTCH_Arch>;
        using strat_ptr = std::shared_ptr<SCOTCH_Strat>;
        using graph_ptr = std::shared_ptr<SCOTCH_Graph>;

    private:
        Integer         m_strategy;
        Matrix          m_A;
        Mat_I           m_node_weights;
        Mat_I           m_edge_weights;
        bool            m_has_weights;
        bool            m_use_edge_weights;

        arch_ptr        m_arch_data;
        graph_ptr       m_graf_data;
        Matrix          m_edge_w;
        Matrix          m_edge_n;
        Integer         m_seed;

    public:
        scotch_wrapper(const Matrix& A, bool use_weights);
        ~scotch_wrapper();

        void            set_strategy(scotch_strategy strategy);
        void            set_seed(Integer s);
        Integer         get_seed() const;

        Matrix          make_assign(Integer n_part, Real lb);
        Matrix          make_assign(Integer n_part, const Matrix& proc_weights, Real lb);
        Matrix          make_coloring();
        void            make_ordering(bool ext, permvec& p, Integer& nblocks, Matrix& blocks,
                            Matrix& tree);
        Matrix          edge_separator(Integer n_part, Real lb);
        Matrix          node_separator(Integer n_part, Real lb);
        Matrix          make_clustering(Integer k, Real lb);

        Matrix          make_assign_fixed(Integer n_part, const Matrix& fixed, Real lb);
        Matrix          make_assign_fixed(Integer n_part, const Matrix& work_shares, 
                            const Matrix& fixed, Real lb);
        Matrix          edge_separator_fixed(Integer n_part, const Matrix& fixed, Real lb);        

        template<class V>
        Matrix          make_assign_impl(Integer n_part, const Matrix& proc_weights, Real lb);
        template<class V>
        Matrix          make_assign_fixed_impl(Integer n_part, const Matrix& proc_weights, 
                            Matrix& fixed, Real lb);

        template<class V>
        Matrix          make_coloring_impl();        
        template<class V>
        void            make_ordering_impl(bool ext, permvec& p, Integer& nblocks, Matrix& blocks,
                            Matrix& tree);

        template<class V>
        Matrix          make_separator_impl(Integer n_part, bool node_sep, bool cluster, Real lb);
        template<class V>
        Matrix          make_separator_fixed_impl(Integer n_part, Matrix& fixed, Real lb);

        void            set_node_weights(const Matrix& W);

    private:
        void            init_arch(std::shared_ptr<SCOTCH_Arch>& arch_data);
        void            init_strategy(std::shared_ptr<SCOTCH_Strat>& strat_data);

        template<class V>
        void            init_graph(std::shared_ptr<SCOTCH_Graph>& graph_data);

        void            check_error(int err, const std::string& msg);        
        Matrix          make_fixed(const Matrix& fixed, Integer n_part);

        void            correct_load(Real& lb) const;
        Integer         default_strategy() const;

        static void     arch_deleter(SCOTCH_Arch*);
        static void     graph_deleter(SCOTCH_Graph*);
        static void     strategy_deleter(SCOTCH_Strat*);

        scotch_wrapper(const scotch_wrapper&) = delete;
        scotch_wrapper& operator=(const scotch_wrapper&) = delete;
};

class scotch_partit_impl : public scotch_wrapper
{
    public:
        scotch_partit_impl(const Matrix& A, bool use_weights)
            :scotch_wrapper(A, use_weights)
        {};
};

}};
