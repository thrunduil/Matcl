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

#include "matcl-linalg/ordering/scotch_wrapper.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/graph/graph_manip.h"
#include "part_wrappers_utils.h"

#pragma warning (push)
#pragma warning (disable:4127)

//#include "petsc/mpiuni/mpi.h"
//#include "petscsys.h"

namespace matcl { namespace details
{

struct scotch_init
{
    scotch_init()
    {
        SCOTCH_installErrorHandler(&error_handler);
    };

    static void error_handler(const char* msg)
    {
        throw error::scotch_error(msg);
    };
};

static scotch_init scotch_initializer;

scotch_graph_wrapper::scotch_graph_wrapper(const Matrix& A, bool use_weights)
    : m_strategy(default_strategy()), m_has_weights(false), m_seed(0)
    , m_use_edge_weights(use_weights), m_A(make_adjancency_matrix(A))
    , m_node_weights(ti::ti_empty()), m_edge_weights(ti::ti_empty())
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    init_arch(m_arch_data);
    SCOTCH_randomSeed(m_seed);

    if (m_use_edge_weights)
        wrappers_utils().build_edge_weights(m_edge_weights, m_A, true, 1000);
};

scotch_graph_wrapper::~scotch_graph_wrapper()
{};

void scotch_graph_wrapper::correct_load(Real& lb) const
{
    if (lb < 0.0)
        lb  = 0.0;
    if (lb > 1.0)
        lb = 1.0;
}

Integer scotch_graph_wrapper::default_strategy() const
{
    return SCOTCH_STRATQUALITY | SCOTCH_STRATSAFETY;    
};

void scotch_graph_wrapper::init_arch(std::shared_ptr<SCOTCH_Arch>& arch_data)
{
    SCOTCH_Arch* arch = new SCOTCH_Arch();
    int err = SCOTCH_archInit(arch);

    if (err)
        delete arch;

    check_error(err, "unable to initialize target architecture");

    arch_data.reset(arch, &arch_deleter);
};

template<class V>
void scotch_graph_wrapper::init_graph(std::shared_ptr<SCOTCH_Graph>& graph_data)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;

    SCOTCH_Graph* graph = new SCOTCH_Graph();

    int err = SCOTCH_graphInit(graph);

    if (err)
        delete graph;

    check_error(err, "unable to initialize graph");    

    graph_data.reset(graph, &graph_deleter);

    const Mat_S& mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    

    int array_base      = 0;    
    int num_edges       = mat_A.nnz();
    Integer* col_start  = mat_A.rep().ptr_c();
    Integer* col_end    = mat_A.rep().ptr_c() + 1;
    Integer* edges      = mat_A.rep().ptr_r();
    int num_nodes       = m_A.rows();

    const int* node_weights = nullptr;
    int* label_array        = nullptr;
    const int* edge_weights = nullptr;

    // build edge weights
    if (m_use_edge_weights == true)
        edge_weights    = m_edge_weights.ptr();

    // build node weights
    if (m_has_weights == true)
        node_weights    = m_node_weights.ptr();

    err = SCOTCH_graphBuild(graph, array_base, num_nodes, col_start, col_end, 
                    node_weights, label_array, num_edges, edges, edge_weights);

    check_error(err, "unable to build graph");

    //SCOTCH_graphCheck(graph);
};

void scotch_graph_wrapper::init_strategy(std::shared_ptr<SCOTCH_Strat>& strat_data)
{
    SCOTCH_Strat* strat = new SCOTCH_Strat();
    int err = SCOTCH_stratInit(strat);

    if (err)
        delete strat;

    check_error(err, "unable to initialize strategy");

    //remove randomness
    SCOTCH_randomReset();

    strat_data.reset(strat, &strategy_deleter);
};

void scotch_graph_wrapper::arch_deleter(SCOTCH_Arch* arch)
{
    if (arch)
    {
        SCOTCH_archExit(arch);
        delete arch;
    };
};

void scotch_graph_wrapper::graph_deleter(SCOTCH_Graph* graph)
{
    if (graph)
    {
        SCOTCH_graphExit(graph);    
        delete graph;
    };
};

void scotch_graph_wrapper::strategy_deleter(SCOTCH_Strat* strat)
{
    if (strat)
    {
        SCOTCH_stratExit(strat);    
        delete strat;
    };
};

void scotch_graph_wrapper::set_node_weights(const Matrix& W)
{
    if (W.is_vector() == false)
        throw error::vector_required(W.rows(), W.cols());

    if (W.length() != m_A.rows())
        throw error::invalid_size2(W.rows(), W.cols(), m_A.rows(), 1);

    wrappers_utils().build_node_weights(m_node_weights, W, true, 1000);
    m_has_weights   = true;
};

// set partitioning strategy
void scotch_graph_wrapper::set_strategy(scotch_strategy strategy)
{
    switch (strategy)
    {
        case scotch_strategy::default_val:
            m_strategy  = SCOTCH_STRATDEFAULT | SCOTCH_STRATSAFETY;
            break;
        case scotch_strategy::balance:
            m_strategy  = SCOTCH_STRATBALANCE | SCOTCH_STRATSAFETY;
            break;
        case scotch_strategy::quality:
            m_strategy  = SCOTCH_STRATQUALITY | SCOTCH_STRATSAFETY;
            break;
        case scotch_strategy::speed:
            m_strategy  = SCOTCH_STRATSPEED | SCOTCH_STRATSAFETY;
            break;
    };
};

void scotch_graph_wrapper::set_seed(Integer s)
{
    m_seed = s;
    SCOTCH_randomSeed(s);
};

Integer scotch_graph_wrapper::get_seed() const
{
    return m_seed;
};

struct visit_make_partition : public md::extract_type_switch<Matrix, visit_make_partition,true>
{
    template<class Mat, class Owner, class ... Args>
    static Matrix eval(const matcl::Matrix&, const Mat&, Owner& owner,
                     Args&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.make_assign_impl<V>(std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static Matrix eval_scalar(const matcl::Matrix&, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.make_assign_impl<Mat>(std::forward<Arg>(args)...);
    };
};

struct visit_make_partition_fixed : public md::extract_type_switch<Matrix, visit_make_partition_fixed,true>
{
    template<class Mat, class Owner, class ... Args>
    static Matrix eval(const matcl::Matrix&, const Mat&, Owner& owner,
                     Args&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.make_assign_fixed_impl<V>(std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static Matrix eval_scalar(const matcl::Matrix&, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.make_assign_fixed_impl<Mat>(std::forward<Arg>(args)...);
    };
};

struct visit_make_coloring : public md::extract_type_switch<Matrix, visit_make_coloring,true>
{
    template<class Mat, class Owner, class ... Args>
    static Matrix eval(const matcl::Matrix&, const Mat&, Owner& owner,
                     Args&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.make_coloring_impl<V>(std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static Matrix eval_scalar(const matcl::Matrix&, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.make_coloring_impl<Mat>(std::forward<Arg>(args)...);
    };
};

struct visit_make_ordering : public md::extract_type_switch<void, visit_make_ordering,true>
{
    template<class Mat, class Owner, class ... Args>
    static void eval(const matcl::Matrix&, const Mat&, Owner& owner,
                     Args&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.make_ordering_impl<V>(std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static void eval_scalar(const matcl::Matrix&, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.make_ordering_impl<Mat>(std::forward<Arg>(args)...);
    };
};

struct visit_make_separator : public md::extract_type_switch<Matrix, visit_make_separator,true>
{
    template<class Mat, class Owner, class ... Args>
    static Matrix eval(const matcl::Matrix&, const Mat&, Owner& owner,
                     Args&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.make_separator_impl<V>(std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static Matrix eval_scalar(const matcl::Matrix&, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.make_separator_impl<Mat>(std::forward<Arg>(args)...);
    };
};

struct visit_make_separator_fixed : public md::extract_type_switch<Matrix, visit_make_separator_fixed,true>
{
    template<class Mat, class Owner, class ... Args>
    static Matrix eval(const matcl::Matrix&, const Mat&, Owner& owner,
                     Args&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.make_separator_fixed_impl<V>(std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static Matrix eval_scalar(const matcl::Matrix&, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.make_separator_fixed_impl<Mat>(std::forward<Arg>(args)...);
    };
};

Matrix scotch_graph_wrapper::make_assign(Integer n_part, Real lb)
{
    Integer N           = m_A.rows();
    n_part              = std::max(1, std::min(n_part, N));
    Matrix proc_weights = repmat(1.0, n_part, 1);

    return make_assign(n_part, proc_weights, lb);
};

Matrix scotch_graph_wrapper::make_assign_fixed(Integer n_part, const Matrix& fixed, Real lb)
{
    Integer N           = m_A.rows();
    n_part              = std::max(1, std::min(n_part, N));
    Matrix proc_weights = repmat(1.0, n_part, 1);

    return make_assign_fixed(n_part, proc_weights, fixed, lb);
};

Matrix scotch_graph_wrapper::make_assign(Integer n_part, const Matrix& proc_weights, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    if (proc_weights.length() != n_part)
        throw error::invalid_size2(proc_weights.rows(), proc_weights.cols(), n_part, 1);
    
    value_code vc0      = m_A.get_value_code();

    if (n_part == 1)
    {
        value_code vc   = matrix_traits::real_value_type(vc0);
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

        return ones(num_nodes, 1, vc);
    };

    Matrix dum = zeros(0,0,vc0);

    return visit_make_partition::make<const Matrix&>(dum, *this, n_part, proc_weights, lb);
};

Matrix scotch_graph_wrapper::make_assign_fixed(Integer n_part, const Matrix& proc_weights, 
                            const Matrix& fixed0, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    if (proc_weights.length() != n_part)
        throw error::invalid_size2(proc_weights.rows(), proc_weights.cols(), n_part, 1);
    
    value_code vc0      = m_A.get_value_code();
    Matrix fixed        = make_fixed(fixed0, n_part);

    if (n_part == 1)
    {
        value_code vc   = matrix_traits::real_value_type(vc0);
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

        return ones(num_nodes, 1, vc);
    };

    Matrix dum = zeros(0,0,vc0);

    return visit_make_partition_fixed::make<const Matrix&>(dum, *this, n_part, proc_weights, fixed, lb);
}
        
Matrix scotch_graph_wrapper::edge_separator(Integer n_part, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();

    if (n_part == 1)
    {
        value_code vc   = matrix_traits::real_value_type(vc0);
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

        return ones(num_nodes, 1, vc);
    };

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator::make<const Matrix&>(dum, *this, n_part, false, false, lb);
};

Matrix scotch_graph_wrapper::make_clustering(Integer n_part, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();

    if (n_part == 1)
    {
        value_code vc   = matrix_traits::real_value_type(vc0);
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

        return ones(num_nodes, 1, vc);
    };

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator::make<const Matrix&>(dum, *this, n_part, false, true, lb);
};

Matrix scotch_graph_wrapper::edge_separator_fixed(Integer n_part, const Matrix& fixed0, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();
    Matrix fixed        = make_fixed(fixed0, n_part);

    if (n_part == 1)
    {
        value_code vc   = matrix_traits::real_value_type(vc0);
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

        return ones(num_nodes, 1, vc);
    };

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator_fixed::make<const Matrix&>(dum, *this, n_part, fixed, lb);
};

Matrix scotch_graph_wrapper::node_separator(Integer n_part, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();

    if (n_part == 1)
    {
        value_code vc   = matrix_traits::real_value_type(vc0);
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

        return ones(num_nodes, 1, vc);
    };

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator::make<const Matrix&>(dum, *this, n_part, true, false, lb);
};

Matrix scotch_graph_wrapper::make_coloring()
{
    value_code vc0      = m_A.get_value_code();
    Matrix dum          = zeros(0,0,vc0);

    return visit_make_coloring::make<const Matrix&>(dum, *this);
};

void scotch_graph_wrapper::make_ordering(bool ext, permvec& p, Integer& nblocks, Matrix& blocks, Matrix& tree)
{
    value_code vc0  = m_A.get_value_code();
    Matrix dum      = zeros(0, 0, vc0);

    return visit_make_ordering::make<const Matrix&>(dum, *this, ext, p, nblocks, blocks, tree);
};

template<class V>
Matrix scotch_graph_wrapper::make_assign_impl(Integer n_part, const Matrix& proc_weights, Real lb)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    const Mat_S& mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    

    int num_proc        = n_part;
    int num_nodes       = m_A.rows();

    Mat_I mat_part_w(ti::ti_empty(), n_part, 1);
    Integer * pw_int    = mat_part_w.ptr();
    const Real* pw_real = proc_weights.get_array<Real>();
    Real pw_sum         = 0.0;

    for (Integer i = 0; i < num_proc; ++i)
        pw_sum          += abs(pw_real[i]);

    pw_sum              = pw_sum / num_proc;
    if (pw_sum == 0.0)
        pw_sum          = 1.0;

    for (Integer i = 0; i < num_proc; ++i)
    {
        // just rescale to obtain value 1 for processor with 0.1% share
        Real val        = std::max(1.0, 1000.0 * (abs(pw_real[i]) / pw_sum));
        pw_int[i]       = (Integer)val;
    }

    Mat_I mat_res(ti::ti_empty(), num_nodes, 1);
    Integer* ptr_res    = mat_res.ptr();    

    if (!m_graf_data)
        init_graph<V>(m_graf_data);

    strat_ptr m_strat_data;
    init_strategy(m_strat_data);

    correct_load(lb);
    int err = SCOTCH_stratGraphMapBuild(m_strat_data.get(), m_strategy, num_proc, lb);
    check_error(err, "unable to build strategy");

    err     = SCOTCH_archCmpltw(m_arch_data.get(), num_proc, pw_int);
    check_error(err, "archCmpltw failed");

    err = SCOTCH_graphMap(m_graf_data.get(), m_arch_data.get(), m_strat_data.get(), ptr_res);

    check_error(err, "unable to compute assigning");

    Matrix ret  = wrappers_utils().convert_to_partitioning(mat_res, num_proc, false, false, 
                                                           m_A.get_value_code());
    return ret;
};

template<class V>
Matrix scotch_graph_wrapper::make_assign_fixed_impl(Integer n_part, const Matrix& proc_weights, 
                            Matrix& fixed, Real lb)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    const Mat_S& mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    

    int num_proc        = n_part;

    Mat_I mat_part_w(ti::ti_empty(), n_part, 1);
    Integer * pw_int    = mat_part_w.ptr();
    const Real* pw_real = proc_weights.get_array<Real>();
    Real pw_sum         = 0.0;

    for (Integer i = 0; i < num_proc; ++i)
        pw_sum          += abs(pw_real[i]);

    pw_sum              = pw_sum / num_proc;
    if (pw_sum == 0.0)
        pw_sum          = 1.0;

    for (Integer i = 0; i < num_proc; ++i)
    {
        // just rescale to obtain value 1 for processor with 0.1% share
        Real val        = std::max(1.0, 1000.0 * (abs(pw_real[i]) / pw_sum));
        pw_int[i]       = (Integer)val;
    }

    Mat_I ifixed        = fixed.impl_unique<Mat_I>();
    Integer* ptr_res    = ifixed.ptr();

    if (!m_graf_data)
        init_graph<V>(m_graf_data);

    strat_ptr m_strat_data;
    init_strategy(m_strat_data);

    int err = SCOTCH_archCmpltw(m_arch_data.get(), num_proc, pw_int);
    check_error(err, "archCmpltw failed");

    correct_load(lb);
    err     = SCOTCH_stratGraphMapBuild(m_strat_data.get(), m_strategy, num_proc, lb);
    check_error(err, "unable to build strategy");

    err = SCOTCH_graphMapFixed(m_graf_data.get(), m_arch_data.get(), m_strat_data.get(), ptr_res);

    check_error(err, "unable to compute assigning");

    Matrix ret  = wrappers_utils().convert_to_partitioning(ifixed, num_proc, false, false,
                                                           m_A.get_value_code());
    return ret;
};

template<class V>
Matrix scotch_graph_wrapper::make_coloring_impl()
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    const Mat_S& mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    
    int num_nodes       = m_A.rows();

    Mat_I mat_res(ti::ti_empty(), num_nodes, 1);
    Integer* ptr_res    = mat_res.ptr();

    if (!m_graf_data)
        init_graph<V>(m_graf_data);

    Integer n_colors;
    int err     = SCOTCH_graphColor(m_graf_data.get(), ptr_res, &n_colors, 0);

    check_error(err, "unable to compute coloring");

    Matrix ret  = wrappers_utils().convert_to_partitioning(mat_res, n_colors, false, true,
                                                           m_A.get_value_code());
    return ret;
};

template<class V>
void scotch_graph_wrapper::make_ordering_impl(bool ext, permvec& p, Integer& nblocks, Matrix& blocks,
                            Matrix& tree)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    const Mat_S& mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    
    int num_nodes       = m_A.rows();

    Mat_I mat_res(ti::ti_empty(), num_nodes, 1);
    Integer* ptr_perm   = nullptr;
    Integer* perm_inv   = mat_res.ptr();
    Integer* ptr_range  = nullptr;
    Integer* ptr_tree   = nullptr;

    if (ext == true)
    {
        blocks          = make_integer_dense_noinit(num_nodes + 1, 1, ptr_range);
        tree            = make_integer_dense_noinit(num_nodes, 1, ptr_tree);
    };

    if (!m_graf_data)
        init_graph<V>(m_graf_data);

    strat_ptr m_strat_ord;
    init_strategy(m_strat_ord);
    
    int num_nested_dis  = 0;    // default in scotch
    Real imbalance      = 0.2;  // default in scotch

    int err     = SCOTCH_stratGraphOrderBuild(m_strat_ord.get(), m_strategy, num_nested_dis, imbalance);
    check_error(err, "unable to build strategy");

    err         = SCOTCH_graphOrder(m_graf_data.get(), m_strat_ord.get(), ptr_perm,
                                    perm_inv, &nblocks, ptr_range, ptr_tree);

    check_error(err, "unable to compute ordering");

    for (Integer i = 0; i < num_nodes; ++i)
        perm_inv[i] += 1;

    p   = permvec::from_matrix(Matrix(mat_res,false));

    if (ext == false)
        return;

    for (Integer i = 0; i < nblocks + 1; ++i)
        ptr_range[i]    += 1;

    for (Integer i = 0; i < nblocks; ++i)
    {
        if (ptr_tree[i] != -1)
            ptr_tree[i] += 1;
    }

    blocks.resize(nblocks + 1, 1);
    tree.resize(nblocks, 1);

    return;
};

template<class V>
Matrix scotch_graph_wrapper::make_separator_impl(Integer n_part, bool node_sep, bool cluster, Real lb)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    const Mat_S& mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    
    int num_nodes       = m_A.rows();

    Mat_I mat_res(ti::ti_empty(), num_nodes, 1);
    Integer* ptr_res    = mat_res.ptr();

    if (!m_graf_data)
        init_graph<V>(m_graf_data);

    strat_ptr m_strat_edge_sep;
    init_strategy(m_strat_edge_sep);

    int err;
 
    int cl_size = 5;
    Real dense  = 0.7;

    correct_load(lb);

    if (node_sep == true)
        err     = SCOTCH_stratGraphPartOvlBuild(m_strat_edge_sep.get(), m_strategy, n_part, lb);
    else if (cluster == true)
        err     = SCOTCH_stratGraphClusterBuild(m_strat_edge_sep.get(), m_strategy, cl_size, dense, lb);
    else
        err     = SCOTCH_stratGraphMapBuild(m_strat_edge_sep.get(), m_strategy, n_part, lb);

    check_error(err, "unable to build strategy");

    if (node_sep == false || cluster == true)
        err    = SCOTCH_graphPart(m_graf_data.get(), n_part, m_strat_edge_sep.get(), ptr_res);
    else
        err    = SCOTCH_graphPartOvl(m_graf_data.get(), n_part, m_strat_edge_sep.get(), ptr_res);

    check_error(err, "unable to compute separator");

    Matrix ret  = wrappers_utils().convert_to_partitioning(mat_res, n_part, node_sep, true,
                                                           m_A.get_value_code());
    return ret;
};

template<class V>
Matrix scotch_graph_wrapper::make_separator_fixed_impl(Integer n_part, Matrix& fixed, Real lb)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V,struct_sparse>;
    const Mat_S& mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    
    Mat_I mat_res       = fixed.impl_unique<Mat_I>();
    Integer* ptr_res    = mat_res.ptr();

    if (!m_graf_data)
        init_graph<V>(m_graf_data);

    strat_ptr m_strat_edge_sep;
    init_strategy(m_strat_edge_sep);

    correct_load(lb);

    int err = SCOTCH_stratGraphMapBuild(m_strat_edge_sep.get(), m_strategy, n_part, lb);
    check_error(err, "unable to build strategy");

    err     = SCOTCH_graphPartFixed(m_graf_data.get(), n_part, m_strat_edge_sep.get(), ptr_res);

    check_error(err, "unable to compute separator");

    Matrix ret  = wrappers_utils().convert_to_partitioning(mat_res, n_part, false, true,
                                                           m_A.get_value_code());
    return ret;

};

void scotch_graph_wrapper::check_error(int err, const std::string& msg)
{
    if (err)
        throw error::scotch_error(msg);
}

Matrix scotch_graph_wrapper::make_fixed(const Matrix& fixed, Integer n_part)
{
    Integer N   = m_A.rows();

    if (fixed.is_vector() == false)
        throw error::vector_required(fixed.rows(), fixed.cols());
    if (fixed.length() != N)
        throw error::invalid_size2(fixed.rows(), fixed.cols(), N, 1);

    Integer* fix_ptr;
    Matrix fix  = make_integer_dense_noinit(N, 1, fix_ptr);

    const Integer* fix_in   = fixed.get_array<Integer>();

    for (Integer i = 0; i < N; ++i)
    {
        Integer val = fix_in[i];

        if (val <= 0)
        {
            fix_ptr[i]  = -1;
        }
        else
        {
            if (val <= n_part)
                fix_ptr[i]  = val - 1;
            else
            {
                std::ostringstream msg;
                msg << "group " << val << " does not exist";
                throw error::scotch_error(msg.str());
            }
        };
    }

    return fix;
};

}};

#pragma warning (pop)