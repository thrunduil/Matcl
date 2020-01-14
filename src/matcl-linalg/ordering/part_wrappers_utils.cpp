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

#include "part_wrappers_utils.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

namespace matcl { namespace details
{

template<class V, bool Is_complex = md::is_complex<V>::value>
struct to_int_rep
{
    static Integer eval(const V& val)
    {
        return (Integer)val;
    }
};

template<class V>
struct to_int_rep<V,true>
{
    static Integer eval(const V& val)
    {
        return (Integer)abs(val);
    }
};

template<>
struct to_int_rep<Object,false>
{
    static Integer eval(const Object& val)
    {
        return cast_integer(val);
    }
};

struct visit_build_weights : public md::extract_type_switch<void, visit_build_weights,true>
{
    using base_type = md::extract_type_switch<void, visit_build_weights,true>;

    template<class Mat, class Owner, class ... Args>
    static void eval(const matcl::Matrix& W, const Mat&, Owner& owner, Args&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.build_node_weights_impl<V>(W, std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static void eval_scalar(const matcl::Matrix& W, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.build_node_weights_impl<Mat>(W, std::forward<Arg>(args)...);
    };

    template<class S, class Owner, class ... Args>
    static void eval(const matcl::Matrix&, const raw::Matrix<Object,S>&, Owner&, Args&& ...)
    {
        throw error::object_value_type_not_allowed("partitioning");
    };

    template<class Owner, class ... Arg>
    static void eval_scalar(const matcl::Matrix&, const Object&, Owner&, Arg&& ...)
    {
        throw error::object_value_type_not_allowed("partitioning");
    };
};
struct visit_edge_build_weights : public md::extract_type_switch<void, visit_edge_build_weights,true>
{
    using base_type = md::extract_type_switch<void, visit_edge_build_weights,true>;

    template<class Mat, class Owner, class ... Args>
    static void eval(const matcl::Matrix& W, const Mat&, Owner&, Args&& ... args)
    {
        using V = typename Mat::value_type;
        return build_edge_weights_impl<V>::eval(W, std::forward<Args>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static void eval_scalar(const matcl::Matrix& W, const Mat&, Owner&, Arg&& ... args)
    {
        return build_edge_weights_impl<Mat>::eval(W, std::forward<Arg>(args)...);
    };

    template<class S, class Owner, class ... Args>
    static void eval(const matcl::Matrix&, const raw::Matrix<Object,S>&, Owner&, Args&& ...)
    {
        throw error::object_value_type_not_allowed("partitioning");
    };

    template<class Owner, class ... Arg>
    static void eval_scalar(const matcl::Matrix&, const Object&, Owner&, Arg&& ...)
    {
        throw error::object_value_type_not_allowed("partitioning");
    };
};

struct visit_get_structure : public md::extract_type_switch<void, visit_get_structure,true>
{
    using base_type = md::extract_type_switch<void, visit_get_structure,true>;

    template<class Mat, class Owner, class ... Arg>
    static void eval(const matcl::Matrix& A, const Mat&, Owner& owner, Arg&& ... args)
    {
        using V = typename Mat::value_type;
        return owner.get_structure_impl<V>(A, std::forward<Arg>(args)...);
    };

    template<class Mat, class Owner, class ... Arg>
    static void eval_scalar(const matcl::Matrix& A, const Mat&, Owner& owner, Arg&& ... args)
    {
        return owner.get_structure_impl<Mat>(A, std::forward<Arg>(args)...);
    };
};

template<class Val>
struct build_edge_weights_impl
{
    using Mat_I = raw::Matrix<Integer, struct_dense>;

    static void eval(const Matrix& W, Mat_I& edge_weights, bool only_pos, Integer max_val)
    {
        using VR            = typename md::real_type<Val>::type;
        using Mat_S         = raw::Matrix<Val,struct_sparse>;

        const Mat_S& rep    = convert(W, Mat_S::matrix_code).get_impl<Mat_S>();

        Integer num_edges   = rep.nnz();

        Mat_I iedge_w(ti::ti_empty(), num_edges, 1);
        edge_weights.assign_to_fresh(iedge_w);

        Integer off         = rep.rep().offset();
        const Val* vals     = rep.rep().ptr_x() + off;

        Integer* ptr_weights= iedge_w.ptr();

        if (vals == nullptr)
            return;

        VR val_max      = abs(vals[0]);
        VR val_min      = val_max;

        for (Integer i = 1; i < num_edges; ++i)
        {
            VR val       = abs(vals[i]);

            if (val < val_min)
                val_min = val;
            if (val > val_max)
                val_max = val;
        };

        if (val_max == VR(0.0))
            val_max     = VR(1.0);

        // assign load = +-max_val to the largest and rescale proportionally
        VR targer_max   = VR(max_val);
        VR scale        = targer_max / val_max;

        //do not take abs; scotch works also with negative weights
        for (Integer i = 0; i < num_edges; ++i)
        {
            Val val         = vals[i] * scale;
            Integer ival    = to_int_rep<Val>::eval(val);
            
            if (ival == 0)
                ival        = 1;

            if (only_pos)
                ptr_weights[i]  = std::max(1,abs(ival));
            else
                ptr_weights[i]  = ival;
        };

        return;
    };

};

template<>
struct build_edge_weights_impl<Integer>
{
    using Mat_I = raw::Matrix<Integer, struct_dense>;

    static void eval(const Matrix& W, Mat_I& edge_weights, bool only_pos, Integer)
    {
        using Val           = Integer;
        using Mat_S         = raw::Matrix<Val,struct_sparse>;

        const Mat_S& rep    = convert(W, Mat_S::matrix_code).get_impl<Mat_S>();

        Integer num_edges   = rep.nnz();

        Mat_I iedge_w(ti::ti_empty(), num_edges, 1);
        edge_weights.assign_to_fresh(iedge_w);

        Integer off         = rep.rep().offset();
        const Val* vals     = rep.rep().ptr_x() + off;

        Integer* ptr_weights    = iedge_w.ptr();

        if (vals == nullptr)
            return;

        if (only_pos)
        {
            for (Integer i = 0; i < num_edges; ++i)
                ptr_weights[i]  = std::max(1,abs(vals[i]));
        }
        else
        {
            for (Integer i = 0; i < num_edges; ++i)
                ptr_weights[i]  = vals[i];
        };
        return;
    };
};

void wrappers_utils::build_node_weights(raw::Matrix<Integer, struct_dense>& node_weights, const Matrix& W,
                                        bool only_nonneg, Integer max_val)
{
    visit_build_weights::make<const Matrix&>(W, *this, node_weights, only_nonneg, max_val);
};

void wrappers_utils::build_edge_weights(Mat_I& edge_weights, const Matrix& A, bool only_pos, Integer max_val)
{
    visit_edge_build_weights::make<const Matrix&>(A, *this, edge_weights, only_pos, max_val);
};
void wrappers_utils::get_structure(const Matrix& A, const Integer*& ptr_c, const Integer*& ptr_r)
{
    visit_get_structure::make<const Matrix&>(A, *this, ptr_c, ptr_r);
};

template<class Val>
void wrappers_utils::get_structure_impl(const Matrix& A, const Integer*& ptr_c, const Integer*& ptr_r)
{
    using Mat_S         = raw::Matrix<Val,struct_sparse>;
    const Mat_S& rep    = convert(A, Mat_S::matrix_code).get_impl<Mat_S>();
    ptr_c               = rep.rep().ptr_c();
    ptr_r               = rep.rep().ptr_r();
};

template<class V>
void wrappers_utils::build_node_weights_impl(const Matrix& W, Mat_I& node_weights, bool only_nonneg,
                                             Integer max_val)
{
    if (std::is_same<V,Integer>::value)
    {
        Matrix Wc       = convert(W, Mat_I::matrix_code); 
        node_weights.assign_to_fresh(Wc.get_impl<Mat_I>());
        return;
    };

    using VR            = typename md::real_type<V>::type;

    Integer M           = W.rows();
    Integer N           = W.cols();

    Mat_I inodes_w(ti::ti_empty(), M, N);
    node_weights.assign_to_fresh(inodes_w);

    Integer* weights_ptr= inodes_w.ptr();
    const VR* ptr_W     = W.get_array<VR>();
    VR val_max          = abs(ptr_W[0]);

    for (Integer i = 1; i < M * N; ++i)
    {
        VR val       = abs(ptr_W[i]);

        if (val > val_max)
            val_max = val;
    };

    // assign load = +-max_val to the largest and rescale proportionally
    VR targer_max   = VR(max_val);
    VR scale        = targer_max / val_max;

    for (Integer i = 0; i < M * N; ++i)
    {
        VR val          = ptr_W[i] * scale;
        Integer ival    = (Integer)val;
            
        if (ival == 0)
            ival        = 1;

        if (only_nonneg)
            weights_ptr[i]  = abs(ival);
        else
            weights_ptr[i]  = ival;
    };

    return;
};

Matrix wrappers_utils::convert_to_partitioning(const Mat_I& mat, Integer n_part, 
                            bool with_sep, bool rem_empty, value_code vc)
{
    vc              = matrix_traits::real_value_type(vc);

    if (vc == value_code::v_float)
        return convert_to_partitioning_impl<Float>(mat, n_part, with_sep, rem_empty);
    else
        return convert_to_partitioning_impl<Real>(mat, n_part, with_sep, rem_empty);
}

template<class V>
Matrix wrappers_utils::convert_to_partitioning_impl(const Mat_I& mat, Integer n_part, 
                            bool with_sep, bool rem_empty)
{
    const Integer* indices  = mat.ptr();
    Integer N               = mat.rows();

    using Mat_S     = raw::Matrix<V, struct_sparse>;
    using ccs_rep   = raw::details::sparse_ccs<V>;
    Integer sep     = with_sep ? 1 : 0;
    Mat_S part(ti::ti_empty(), N, n_part + sep, N);

    ccs_rep rep     = part.rep();    

    Integer* d_c    = rep.ptr_c();
    Integer* d_r    = rep.ptr_r();
    V* d_x          = rep.ptr_x();

    //count nz
    for (Integer i = 0; i < N; ++i)
        d_c[indices[i] + sep] += 1;

    //create column indices
    Integer nz      = 0;
    for (Integer i = 0; i < n_part + sep; ++i)
    {
        Integer loc_nz  = d_c[i];
        d_c[i]          = nz;
        nz              += loc_nz;
    };
    d_c[n_part + sep]   = nz;
        
    //fill arrays
    for (Integer i = 0; i < N; ++i)
    {
        Integer ind     = indices[i] + sep;
        Integer pos     = d_c[ind];
        d_r[pos]        = i;
        d_x[pos]        = 1.0;
        d_c[ind]        += 1;
    }

    //restore column indices
    for (Integer i = n_part + sep - 1; i >= 1; --i)
        d_c[i]          = d_c[i-1];

    d_c[0]              = 0;

    //remove empty groups, but keep empty separator
    Integer K           = sep;

    if (rem_empty == true)
    {
        nz              = d_c[K+1];
        for (Integer i = K+1; i <= n_part + sep; ++i)
        {
            nz          = d_c[i];
            bool empty  = (d_c[i] == d_c[i-1]);
            d_c[K+1]    = nz;

            if (empty == true)
                continue;

            ++K;
        };
        d_c[K]  = d_c[n_part + sep];
    }
    else
    {
        K   = n_part + sep;
    };

    Matrix ret = Matrix(part,false);

    if (K < n_part + sep)
        ret = ret(colon(), colon(1,K));

    return ret;
};

}};

#pragma warning(pop)