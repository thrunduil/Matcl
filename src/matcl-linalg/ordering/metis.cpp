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

#include "matcl-linalg/graph/metis.h"
#include "matcl-linalg/ordering/scotch_wrapper.h"
#include "matcl-linalg/graph/graph_manip.h"
#include "part_wrappers_utils.h"

#pragma warning (push)
#pragma warning (disable: 4005)
#include "extern/metis/include/metis.h"
#pragma warning (pop)

namespace matcl { namespace details
{

static const Integer MAX_VAL    = 1000;

class metis_impl
{
    private:
        using Mat_I         = raw::Matrix<Integer, struct_dense>;
        using metis_real    = real_t;

    private:
        Matrix              m_A;
        Integer             m_seed;
        const Integer*      m_ptr_c;
        const Integer*      m_ptr_r;
        Mat_I               m_options;
        Mat_I               m_node_weights;
        Mat_I               m_node_size;
        Mat_I               m_edge_weights;
        bool                m_use_edge_weights;
        bool                m_has_weights;
        bool                m_has_size;
        Integer             m_objective;

    public:
        metis_impl(const Matrix& A, bool use_weights);
        ~metis_impl();

        void                set_seed(Integer s);
        Integer             get_seed() const;
        void                set_node_weights(const Matrix& W);
        void                set_node_size(const Matrix& W);
        permvec             ordering();
        Matrix              partition(Integer k, part_type type);
        Matrix              partition(Integer k, part_type type, const Matrix& contraints, 
                                const Matrix& balance);
        Integer*            get_option_ptr();
        const Integer*      get_option_ptr() const;
        Integer             get_objective() const;

    private:
        void                initialize();
        const Integer*      get_column_indices() const;
        const Integer*      get_row_indices() const;
        const Integer*      get_edge_weights() const;
        const Integer*      get_node_weights() const;
        const Integer*      get_node_size() const;
        Integer*            get_options();

        Matrix              make_partition_matrix(const Mat_I& partition, Integer k) const;
};

metis_impl::metis_impl(const Matrix& A, bool use_weights)
    : m_seed(0), m_A(make_adjancency_matrix(A)), m_options(ti::ti_empty())
    , m_node_weights(ti::ti_empty()), m_edge_weights(ti::ti_empty()), m_node_size(ti::ti_empty())
    , m_use_edge_weights(use_weights), m_has_weights(false), m_has_size(false)
    , m_objective(-1)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    initialize();
};

metis_impl::~metis_impl()
{};

Integer metis_impl::get_objective() const
{
    return m_objective;
};

permvec metis_impl::ordering()
{
    Integer N               = m_A.rows();
    const Integer* ptr_c    = get_column_indices();
    const Integer* ptr_r    = get_row_indices();
    const Integer* ptr_w    = get_node_weights();
    const Integer* opt      = get_options();

    Mat_I perm(ti::ti_empty(), N, 1);
    Mat_I iperm(ti::ti_empty(), N, 1);

    Integer* ptr_perm       = perm.ptr();
    Integer* ptr_iperm      = iperm.ptr();

    int err = METIS_NodeND(&N, (idx_t*)ptr_c, (idx_t*)ptr_r, (idx_t*)ptr_w, (idx_t*)opt, 
                           ptr_perm, ptr_iperm);

    if (err != METIS_OK)
        throw error::metis_error(err);

    for (Integer i = 0; i < N; ++i)
        ptr_perm[i]         += 1;

    permvec p   = permvec::from_matrix(Matrix(perm,false));
    return p;
};

Matrix metis_impl::partition(Integer k, part_type type)
{
    m_objective                 = -1;
    Integer N                   = m_A.rows();
    Integer n_constr            = m_has_weights? m_node_weights.rows() : 1;
    const Integer* vsize        = get_node_size();
    const Integer* ptr_c        = get_column_indices();
    const Integer* ptr_r        = get_row_indices();
    const Integer* ptr_nw       = get_node_weights();
    const Integer* ptr_ew       = get_edge_weights();
    const metis_real* ptr_pw    = nullptr;
    const metis_real* ptr_load  = nullptr;
    Integer* opt                = get_options();

    switch (type)
    {
        case part_type::edge_cut_recursive:        
        case part_type::edge_cut_k_way:
            opt[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
            break;        
        case part_type::comm_volume_k_way:
            opt[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_VOL;
            break;
    };

    Mat_I partition(ti::ti_empty(), N, 1);

    Integer* ptr_part           = partition.ptr();

    int err = 0;

    switch (type)
    {
        case part_type::edge_cut_recursive:
            err = METIS_PartGraphRecursive(&N, &n_constr, (idx_t*)ptr_c, (idx_t*)ptr_r, 
                    (idx_t*)ptr_nw, (idx_t*)vsize, (idx_t*)ptr_ew, &k, (metis_real*)ptr_pw,
                    (metis_real*)ptr_load, (idx_t*)opt, &m_objective, ptr_part);
            break;
        case part_type::edge_cut_k_way:
        case part_type::comm_volume_k_way:
            err = METIS_PartGraphKway(&N, &n_constr, (idx_t*)ptr_c, (idx_t*)ptr_r, 
                    (idx_t*)ptr_nw, (idx_t*)vsize, (idx_t*)ptr_ew, &k, (metis_real*)ptr_pw, 
                    (metis_real*)ptr_load, (idx_t*)opt, &m_objective, ptr_part);
            break;
    };

    if (err != METIS_OK)
        throw error::metis_error(err);

    return make_partition_matrix(partition, k);
};

Matrix metis_impl::partition(Integer k, part_type type, const Matrix& contraints, 
                             const Matrix& balance)
{
    m_objective                 = -1;
    Integer N                   = m_A.rows();
    Integer n_constr            = m_has_weights? m_node_weights.rows() : 1;

    if (balance.is_vector() == false)
        throw error::vector_required(balance.rows(), balance.cols());
    if (balance.length() != n_constr)
        throw error::invalid_size2(balance.rows(), balance.cols(), m_A.rows(), 1);

    if (contraints.rows() != n_constr || contraints.cols() != k)
        throw error::invalid_size2(contraints.rows(), contraints.cols(), n_constr, k);
    
    Matrix constr               = abs(contraints) + constants::eps<Float>();
    Matrix sum                  = matcl::sum(constr, 2);
    constr                      = scale_rows(constr, div(1.0, sum));

    using Mat_F = raw::Matrix<Float,struct_dense>;
    Mat_F constr_f(ti::ti_empty(), constr.rows(), constr.cols());
    Mat_F balance_f(ti::ti_empty(), n_constr, 1);

    Float* ptr_constr           = constr_f.ptr();
    Float* ptr_balance          = balance_f.ptr();
    const Real* ptr_constr_mat  = constr.get_array<Real>();
    const Real* ptr_balance_mat = balance.get_array<Real>();

    for (Integer i = 0; i < constr_f.size(); ++i)
        ptr_constr[i]           = (Float)ptr_constr_mat[i];

    for (Integer i = 0; i < n_constr; ++i)
        ptr_balance[i]          = std::max(Float(0.0), (Float)ptr_balance_mat[i]) + 1.0f;

    const Integer* vsize        = get_node_size();
    const Integer* ptr_c        = get_column_indices();
    const Integer* ptr_r        = get_row_indices();
    const Integer* ptr_nw       = get_node_weights();
    const Integer* ptr_ew       = get_edge_weights();
    const metis_real* ptr_pw    = ptr_constr;
    const metis_real* ptr_load  = ptr_balance;
    Integer* opt                = get_options();

    switch (type)
    {
        case part_type::edge_cut_recursive:        
        case part_type::edge_cut_k_way:
            opt[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
            break;        
        case part_type::comm_volume_k_way:
            opt[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_VOL;
            break;
    };

    Mat_I partition(ti::ti_empty(), N, 1);

    Integer* ptr_part           = partition.ptr();

    int err = 0;

    switch (type)
    {
        case part_type::edge_cut_recursive:
            err = METIS_PartGraphRecursive(&N, &n_constr, (idx_t*)ptr_c, (idx_t*)ptr_r, 
                    (idx_t*)ptr_nw, (idx_t*)vsize, (idx_t*)ptr_ew, &k, (metis_real*)ptr_pw,
                    (metis_real*)ptr_load, (idx_t*)opt, &m_objective, ptr_part);
            break;
        case part_type::edge_cut_k_way:
        case part_type::comm_volume_k_way:
            err = METIS_PartGraphKway(&N, &n_constr, (idx_t*)ptr_c, (idx_t*)ptr_r, 
                    (idx_t*)ptr_nw, (idx_t*)vsize, (idx_t*)ptr_ew, &k, (metis_real*)ptr_pw, 
                    (metis_real*)ptr_load, (idx_t*)opt, &m_objective, ptr_part);
            break;
    };

    if (err != METIS_OK)
        throw error::metis_error(err);

    return make_partition_matrix(partition, k);
};

const Integer* metis_impl::get_column_indices() const
{
    return m_ptr_c;
};
const Integer* metis_impl::get_row_indices() const
{
    return m_ptr_r;
};
const Integer* metis_impl::get_edge_weights() const
{
    return m_use_edge_weights? m_edge_weights.ptr() : nullptr;
};
const Integer* metis_impl::get_node_weights() const
{
    return m_has_weights? m_node_weights.ptr() : nullptr;
};
const Integer* metis_impl::get_node_size() const
{
    return m_has_size? m_node_size.ptr() : nullptr;
};

Integer* metis_impl::get_options()
{
    Integer* opts           = m_options.ptr();
    opts[METIS_OPTION_SEED] = m_seed;

    return opts;
};

void metis_impl::initialize()
{
    if (m_use_edge_weights == true)
    {
        wrappers_utils().build_edge_weights(m_edge_weights, m_A, true, MAX_VAL);
    };

    wrappers_utils().get_structure(m_A, m_ptr_c, m_ptr_r);

    //init options
    m_options.assign_to_fresh(Mat_I(ti::ti_empty(), METIS_NOPTIONS, 1));
    Integer* opts   = m_options.ptr();

    METIS_SetDefaultOptions(opts);

    opts[METIS_OPTION_DBGLVL]   = 0;
    opts[METIS_OPTION_SEED]     = 0;
    opts[METIS_OPTION_NUMBERING]= 0;
};

void metis_impl::set_seed(Integer s)
{
    m_seed = s;
}
Integer metis_impl::get_seed() const
{
    return m_seed;
};

Integer* metis_impl::get_option_ptr()
{
    return m_options.ptr();
};
const Integer* metis_impl::get_option_ptr() const
{
    return m_options.ptr();
};

void metis_impl::set_node_weights(const Matrix& W)
{
    Integer K   = W.rows();
    Integer N   = W.cols();

    if (N != m_A.rows())
        throw error::invalid_size2(W.rows(), W.cols(), W.rows(), m_A.rows());

    m_has_weights   = false;

    if (K == 0)
        return;
        
    wrappers_utils().build_node_weights(m_node_weights, W, true, MAX_VAL);
    m_has_weights   = true;
};

void metis_impl::set_node_size(const Matrix& W)
{
    if (W.is_vector() == false)
        throw error::vector_required(W.rows(), W.cols());

    if (W.length() != m_A.rows())
        throw error::invalid_size2(W.rows(), W.cols(), m_A.rows(), 1);

    m_has_size  = false;

    wrappers_utils().build_node_weights(m_node_size, W, true, MAX_VAL);
    m_has_weights   = true;
};

Matrix metis_impl::make_partition_matrix(const Mat_I& partition, Integer k) const
{
    return wrappers_utils().convert_to_partitioning(partition, k, false, false, m_A.get_value_code());
};

}};

namespace matcl
{

metis::metis()
    :m_impl(new details::metis_impl(1.0, false))
{}
metis::metis(const Matrix& A, bool use_weights)
    :m_impl(new details::metis_impl(A,use_weights))
{}

metis::~metis()
{};

permvec metis::ordering()
{
    return m_impl->ordering();
};

Matrix metis::partition(Integer k, part_type type)
{
    return m_impl->partition(k, type);
};
Matrix metis::partition(Integer k, const Matrix& contraints, const Matrix& balance,
                    part_type type)
{
    return m_impl->partition(k, type, contraints, balance);
};
void metis::set_node_weights(const Matrix& W)
{
    return m_impl->set_node_weights(W);
};

void metis::set_node_size(const Matrix& W)
{
    return m_impl->set_node_size(W);
};

void metis::set_seed(Integer s)
{
    return m_impl->set_seed(s);
};

Integer metis::get_seed() const
{
    return m_impl->get_seed();
};

Integer* metis::get_option_ptr()
{
    return m_impl->get_option_ptr();
};
const Integer* metis::get_option_ptr() const
{
    return m_impl->get_option_ptr();
};

Integer metis::objective() const
{
    return m_impl->get_objective();
};


metis& metis::set_coarsing_alg(coarsing_type type)
{
    switch(type)
    {
        case coarsing_type::random:
            m_impl->get_option_ptr()[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
            break;
        case coarsing_type::hem:
            m_impl->get_option_ptr()[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
            break;
    };
    return *this;
};
metis& metis::set_initial_partitioning(initial_partitioning type)
{
    switch(type)
    {
        case initial_partitioning::greedy:
            m_impl->get_option_ptr()[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
            break;
        case initial_partitioning::random:
            m_impl->get_option_ptr()[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM ;
            break;
        case initial_partitioning::edge:
            m_impl->get_option_ptr()[METIS_OPTION_IPTYPE] = METIS_IPTYPE_EDGE;
            break;
        case initial_partitioning::node:
            m_impl->get_option_ptr()[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
            break;
    };
    return *this;
};

metis& metis::set_refinement_type(refinement_type type)
{
    switch(type)
    {
        case refinement_type::fm:
            m_impl->get_option_ptr()[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
            break;
        case refinement_type::greedy:
            m_impl->get_option_ptr()[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY ;
            break;
        case refinement_type::fm_2_sided:
            m_impl->get_option_ptr()[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED;
            break;
        case refinement_type::fm_1_sided:
            m_impl->get_option_ptr()[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP1SIDED;
            break;
    };
    return *this;
};

metis& metis::set_n_cuts(Integer n)
{
    m_impl->get_option_ptr()[METIS_OPTION_NCUTS] = std::max(1, n);
    return *this;
}
metis& metis::set_n_sep(Integer n)
{
    m_impl->get_option_ptr()[METIS_OPTION_NSEPS] = std::max(1, n);
    return *this;
}
metis& metis::set_iter(Integer n)
{
    m_impl->get_option_ptr()[METIS_OPTION_NITER] = std::max(0, n);
    return *this;
}
metis& metis::set_min_conn(bool val)
{
    m_impl->get_option_ptr()[METIS_OPTION_MINCONN] = val ? 1 : 0;
    return *this;
}
metis& metis::set_no_2hop(bool val)
{
    m_impl->get_option_ptr()[METIS_OPTION_NO2HOP] = val ? 1 : 0;
    return *this;
}
metis& metis::set_contig(bool val)
{
    m_impl->get_option_ptr()[METIS_OPTION_CONTIG] = val ? 1 : 0;
    return *this;
}
metis& metis::set_compress(bool val)
{
    m_impl->get_option_ptr()[METIS_OPTION_COMPRESS] = val ? 1 : 0;
    return *this;
}
metis& metis::set_identify_connected(bool val)
{
    m_impl->get_option_ptr()[METIS_OPTION_CCORDER] = val ? 1 : 0;
    return *this;
}
metis& metis::set_pfactor(Real val)
{
    Integer ival    = std::max(Integer(val * 10), 0);
    m_impl->get_option_ptr()[METIS_OPTION_PFACTOR] = ival;
    return *this;
}
metis& metis::set_imbalance(Real val)
{
    Integer ival    = std::max(Integer(val * 1000), 0);
    m_impl->get_option_ptr()[METIS_OPTION_UFACTOR] = ival;
    return *this;
}

};