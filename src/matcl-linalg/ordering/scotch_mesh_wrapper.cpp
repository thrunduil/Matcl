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

namespace matcl { namespace details
{

scotch_mesh_wrapper::scotch_mesh_wrapper(const Matrix& A)
    : m_A(A)
{
    SCOTCH_randomSeed(m_seed);
};

scotch_mesh_wrapper::~scotch_mesh_wrapper()
{};

void scotch_mesh_wrapper::set_seed(Integer s)
{
    m_seed = s;
    SCOTCH_randomSeed(s);
};

Integer scotch_mesh_wrapper::get_seed() const
{
    return m_seed;
};

struct visit_eval_mat : public md::extract_type_switch<void, visit_eval_mat, true>
{
    template<class Mat, class Evaler, class ... Args>
    static void eval(const matcl::Matrix&, const Mat& m, const Evaler& fun, Args&& ... args)
    {
        return fun(m, std::forward<Args>(args)...);
    };

    template<class Scal, class Evaler, class ... Args>
    static void eval_scalar(const matcl::Matrix&, const Scal& s, const Evaler& fun, Args&& ... args)
    {
        using Mat   = raw::Matrix<Scal, struct_dense>;
        Mat Ac      = raw::converter<Mat, Scal>::eval(s);

        return fun(Ac, std::forward<Args>(args)...);
    };
};

template<class Evaler>
void eval_mat(const Matrix& m, const Evaler& fun)
{    
    return visit_eval_mat::make<const Matrix&>(m, fun);
};

void scotch_mesh_wrapper::make_ordering()
{
    value_code vc0  = m_A.get_value_code();
    Matrix dum      = zeros(0,0,vc0);

    auto fun        = [this](auto mat) {return this->make_ordering_impl(mat); };

    return eval_mat(dum, fun);
};

void scotch_mesh_wrapper::check_error(int err, const std::string& msg)
{
    if (err)
        throw error::scotch_error(msg);
}

void scotch_mesh_wrapper::mesh_deleter(SCOTCH_Mesh* mesh)
{
    if (mesh)
    {
        SCOTCH_meshExit(mesh);    
        delete mesh;
    };
};

template<class V>
void scotch_mesh_wrapper::init_mesh(std::shared_ptr<SCOTCH_Mesh>& mesh_data)
{
    using VR            = typename md::real_type<V>::type;
    using Mat_S         = raw::Matrix<V, struct_sparse>;

    SCOTCH_Mesh* mesh   = new SCOTCH_Mesh();
    int err             = SCOTCH_meshInit(mesh);

    if (err)
        delete mesh;

    check_error(err, "unable to initialize mesh");    

    mesh_data.reset(mesh, &mesh_deleter);

    Mat_S mat_A0 = m_A.impl<Mat_S>();

    // scotch only reads this matrix; we can remove const
    Mat_S mat_A         = mat_A0;    

    // min(elem_base, node_base) = 0
    // elem_base = node_num - elem_num or node_base = elem_num - node_num
    int elem_num        = m_A.cols();   // velmnbr
    int node_num        = m_A.rows();   // vnodnbr
    int elem_base       = 0;            // velmbas
    int node_base       = 0;            // vnodbas

    if (node_num >= elem_num)
    {
        elem_base       = node_num - elem_num;
        node_base       = 0;
    }
    else
    {
        node_base       = elem_num - node_num;
        elem_base       = 0;
    };    

    err = SCOTCH_meshBuild(mesh, elem_base, node_base, elem_num, node_num);

    check_error(err, "unable to build graph");

    /*
int SCOTCH meshBuild (,
const SCOTCH Num ,
const SCOTCH Num ,
const SCOTCH Num ,
const SCOTCH Num ,
const SCOTCH Num * verttab,
const SCOTCH Num * vendtab,
const SCOTCH Num * velotab,
const SCOTCH Num * vnlotab,
const SCOTCH Num * vlbltab,
const SCOTCH Num edgenbr,
const SCOTCH Num * edgetab)
*/
    
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

    //SCOTCH_graphCheck(graph);
};

template<class Mat>
void scotch_mesh_wrapper::make_ordering_impl(const Mat&)
{
    using V             = typename Mat::value_type;

    if (!m_mesh_data)
        init_mesh<V>(m_mesh_data);

    //TODO
    /*    
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
    */
};

//TODO
#if 0

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
    value_code vc       = matrix_traits::real_value_type(vc0);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (n_part == 1)
        return ones(num_nodes, 1, vc);

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
    value_code vc       = matrix_traits::real_value_type(vc0);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    Matrix fixed        = make_fixed(fixed0, n_part);

    if (n_part == 1)
        return ones(num_nodes, 1, vc);

    Matrix dum = zeros(0,0,vc0);

    return visit_make_partition_fixed::make<const Matrix&>(dum, *this, n_part, proc_weights, fixed, lb);
}
        
Matrix scotch_graph_wrapper::edge_separator(Integer n_part, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();
    value_code vc       = matrix_traits::real_value_type(vc0);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (n_part == 1)
        return ones(num_nodes, 1, vc);

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator::make<const Matrix&>(dum, *this, n_part, false, false, lb);
};

Matrix scotch_graph_wrapper::make_clustering(Integer n_part, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();
    value_code vc       = matrix_traits::real_value_type(vc0);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (n_part == 1)
        return ones(num_nodes, 1, vc);

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator::make<const Matrix&>(dum, *this, n_part, false, true, lb);
};

Matrix scotch_graph_wrapper::edge_separator_fixed(Integer n_part, const Matrix& fixed0, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();
    value_code vc       = matrix_traits::real_value_type(vc0);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    Matrix fixed        = make_fixed(fixed0, n_part);

    if (n_part == 1)
        return ones(num_nodes, 1, vc);

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator_fixed::make<const Matrix&>(dum, *this, n_part, fixed, lb);
};

Matrix scotch_graph_wrapper::node_separator(Integer n_part, Real lb)
{
    int num_nodes       = m_A.rows();
    n_part              = std::min(std::max(1, n_part), num_nodes);

    value_code vc0      = m_A.get_value_code();
    value_code vc       = matrix_traits::real_value_type(vc0);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (n_part == 1)
        return ones(num_nodes, 1, vc);

    Matrix dum = zeros(0,0,vc0);

    return visit_make_separator::make<const Matrix&>(dum, *this, n_part, true, false, lb);
};

Matrix scotch_graph_wrapper::make_coloring()
{
    value_code vc0      = m_A.get_value_code();
    value_code vc       = matrix_traits::real_value_type(vc0);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    Matrix dum = zeros(0,0,vc0);

    return visit_make_coloring::make<const Matrix&>(dum, *this);
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

#endif

}};

#pragma warning (pop)