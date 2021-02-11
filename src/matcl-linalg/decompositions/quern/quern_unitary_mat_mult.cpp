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

#include "matcl-linalg/decompositions/quern/quern_unitary_mat_mult.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/decompositions/quern/include/quern.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/base/utils.h"
#include "matcl-linalg/decompositions/qr_utils.h"

namespace matcl { namespace details
{

namespace mdyp = matcl::dynamic::functions;

template<class Val, class VM, class Struct>
struct quern_mult_struct {};

//---------------------------------------------------------------
//                      DENSE
//---------------------------------------------------------------
template<class Val, class VM>
struct quern_mult_struct<Val,VM,struct_dense>
{
    using Mat           = raw::Matrix<VM,struct_dense>;
    using unitary_data  = quern_unitary_mat_mult<Val>;

    static const bool is_real_umatrix   = md::is_float_real_scalar<Val>::value;

    static void eval_right(Matrix& ret, const Mat& X, const unitary_data& ud,
                           trans_type t_unit)
    {
        using VTR       = typename md::unify_types<Val,VM>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;

        Integer crs_rows    = ud.m_rows;

        Integer M, K;

        if (t_unit == trans_type::no_trans)
        {
            M   = ud.m_rows;
            K   = ud.m_cols;
        }
        else
        {
            M   = ud.m_cols;
            K   = ud.m_rows;
        };

        Integer N   = X.cols();

        using VR = typename md::real_type<Val>::type;

        ti::ti_type<Val> tiu;
        ti::ti_type<VTR> ret_ti = ti::get_return_ti<ti::ti_type<VTR>>(mdyp::op_mul::eval(),
                                                                        tiu, ti::get_ti(X));

        if (t_unit == trans_type::no_trans)
        {
            ret_type Y      = make_copy<VTR, VM>::eval(X, crs_rows, X.cols());
            VTR* ptr_Y      = Y.ptr();
            Integer Y_ld    = Y.ld();

            VR  c;
            Val s;

            for (Integer xc = 0; xc < N; ++xc)
            {
                for (Integer i = crs_rows - 1; i >= 0; --i)
                {
                    Integer j_max   = ud.row_start[i+1] - 1;
                    Integer j_min   = ud.row_start[i];

                    for (Integer j = j_max; j >= j_min; --j)
                    {
                        Integer k   = ud.column_index[j];
                        quern::make_val_givens<Val>::decode(ud.value, j, c, s);

                        Val sc = conj(s);

                        if (quern::make_val_givens<Val>::is_swap(c) == true)
                        { 
                            std::swap(ptr_Y[i], ptr_Y[k]);
                        }
                        else
                        {
                            VTR newxk   = c*ptr_Y[k] - s*ptr_Y[i];
                            ptr_Y[i]    = c*ptr_Y[i] + sc*ptr_Y[k];
                            ptr_Y[k]    = newxk;
                        };
                    };
                };

                ptr_Y += Y_ld;
            };

            struct_flag su;
            if (M == K)
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unit,trans_type::no_trans,
                                is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));

            ret = Matrix(Y,true);
            return;
        }
        else
        {
            bool conj_U     = (t_unit == trans_type::trans);
            ret_type Y      = make_copy<VTR, VM>::eval(X, crs_rows, X.cols());
            VTR* ptr_Y      = Y.ptr();
            Integer Y_ld    = Y.ld();

            VR  c;
            Val s;
    
            for (Integer xc = 0; xc < N; ++xc)
            {
                for (Integer i = 0; i < crs_rows; ++i)
                {
                    Integer j_start = ud.row_start[i]; 
                    Integer j_end   = ud.row_start[i+1];

                    for (Integer j = j_start; j < j_end; ++j)
                    {
                        Integer k   = ud.column_index[j];
                        quern::make_val_givens<Val>::decode(ud.value,j, c, s);

                        Val sc = conj(s);

                        if (conj_U)
                            std::swap(sc,s);

                        if (c == VR(1))
                        { 
                            std::swap(ptr_Y[i], ptr_Y[k]);
                        }
                        else
                        {
                            VTR newxk   = c*ptr_Y[k] + s*ptr_Y[i];
                            ptr_Y[i]    = c*ptr_Y[i] - sc*ptr_Y[k];
                            ptr_Y[k]    = newxk;
                        }
                    }
                }

                ptr_Y   += Y_ld;
            };

            Y.assign_to_fresh(Y.resize(M,N));

            struct_flag su;
            if (M == K)
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unit,trans_type::no_trans,
                                is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));

            ret = Matrix(Y,true);

            return;
        };
    };

    // Y = X * Q => Y^T = Q^T * X^T
    static void eval_left(Matrix& ret, const Mat& X, const unitary_data& ud,trans_type t_unit)
    {
        using VTR       = typename md::unify_types<Val,VM>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;

        Integer crs_rows    = ud.m_rows;

        Integer M           = X.rows();
        Integer N           = (t_unit == trans_type::no_trans) ? ud.m_cols : ud.m_rows;

        using VR = typename md::real_type<Val>::type;

        ti::ti_type<Val> tiu;
        ti::ti_type<VTR> ret_ti = ti::get_return_ti<ti::ti_type<VTR>>(mdyp::op_mul::eval(),
                                                                        tiu, ti::get_ti(X));

        if (t_unit == trans_type::no_trans)
        {
            ret_type Y      = make_copy<VTR, VM>::eval(X, M, X.cols());
            VTR* ptr_Y      = Y.ptr();
            Integer Y_ld    = Y.ld();

            VR  c;
            Val s;
    
            for (Integer i = 0; i < crs_rows; ++i)
            {
                Integer j_start = ud.row_start[i]; 
                Integer j_end   = ud.row_start[i+1];

                for (Integer j = j_start; j < j_end; ++j)
                {
                    Integer k   = ud.column_index[j];
                    quern::make_val_givens<Val>::decode(ud.value,j, c, s);

                    VTR* Y_i        = ptr_Y + i * Y_ld;
                    VTR* Y_k        = ptr_Y + k * Y_ld;

                    if (c == VR(1))
                    { 
                        for (Integer l = 0; l < M; ++l)
                            std::swap(Y_i[l], Y_k[l]);
                    }
                    else
                    {
                        Val sc = conj(s);
                        for (Integer l = 0; l < M; ++l)
                        {
                            VTR newxk   = c*Y_k[l] + sc*Y_i[l];
                            Y_i[l]      = c*Y_i[l] - s*Y_k[l];
                            Y_k[l]      = newxk;
                        };
                    }
                }
            }

            Y.assign_to_fresh(Y.resize(M,N));

            struct_flag su;
            if (ud.m_rows == ud.m_cols)
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),su,trans_type::no_trans,t_unit,
                                is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

            ret = Matrix(Y,true);

            return;
        }
        else
        {
            bool conj_U     = (t_unit == trans_type::conj_trans);
            ret_type Y      = make_copy<VTR, VM>::eval(X, M, N);
            VTR* ptr_Y      = Y.ptr();
            Integer Y_ld    = Y.ld();

            VR  c;
            Val s;

            for (Integer i = crs_rows - 1; i >= 0; --i)
            {
                Integer j_max   = ud.row_start[i+1] - 1;
                Integer j_min   = ud.row_start[i];

                for (Integer j = j_max; j >= j_min; --j)
                {
                    Integer k   = ud.column_index[j];
                    quern::make_val_givens<Val>::decode(ud.value, j, c, s);

                    VTR* Y_i        = ptr_Y + i * Y_ld;
                    VTR* Y_k        = ptr_Y + k * Y_ld;

                    if (quern::make_val_givens<Val>::is_swap(c) == true)
                    { 
                        for (Integer l = 0; l < M; ++l)
                            std::swap(Y_i[l], Y_k[l]);
                    }
                    else
                    {
                        Val sc          = conj(s);

                        if (conj_U == true)
                            std::swap(s,sc);

                        for (Integer l = 0; l < M; ++l)
                        {
                            VTR newxk   = c*Y_k[l] - s*Y_i[l];
                            Y_i[l]      = c*Y_i[l] + sc*Y_k[l];
                            Y_k[l]      = newxk;
                        };
                    };
                };
            };

            struct_flag su;
            if (ud.m_cols == ud.m_rows)
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),su,trans_type::no_trans,t_unit,
                                is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

            ret = Matrix(Y,true);
            return;
        }
    };
};

//---------------------------------------------------------------
//                      SPARSE
//---------------------------------------------------------------
template<class Val, class VM>
struct quern_mult_struct<Val,VM,struct_sparse>
{
    using Mat           = raw::Matrix<VM,struct_sparse>;
    using unitary_data  = quern_unitary_mat_mult<Val>;

    static void eval_right(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        using DM    = raw::Matrix<VM,struct_dense>;
        DM full_mat = raw::converter<DM,Mat>::eval(mat);
        return quern_mult_struct<Val, VM, struct_dense>::eval_right(ret,full_mat,ud,t_unit);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        using DM    = raw::Matrix<VM,struct_dense>;
        DM full_mat = raw::converter<DM,Mat>::eval(mat);
        return quern_mult_struct<Val, VM, struct_dense>::eval_left(ret,full_mat,ud,t_unit);
    };
};

//---------------------------------------------------------------
//                      BAND
//---------------------------------------------------------------
template<class Val, class VM>
struct quern_mult_struct<Val,VM,struct_banded>
{
    using Mat           = raw::Matrix<VM,struct_banded>;
    using unitary_data  = quern_unitary_mat_mult<Val>;

    static void eval_right(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        using DM    = raw::Matrix<VM,struct_dense>;
        DM full_mat = raw::converter<DM,Mat>::eval(mat);
        return quern_mult_struct<Val, VM, struct_dense>::eval_right(ret,full_mat,ud,t_unit);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        using DM    = raw::Matrix<VM,struct_dense>;
        DM full_mat = raw::converter<DM,Mat>::eval(mat);
        return quern_mult_struct<Val, VM, struct_dense>::eval_left(ret,full_mat,ud,t_unit);
    };
};

//---------------------------------------------------------------
//              quern_mult_impl
//---------------------------------------------------------------
template<class Val, class VM, class Struct>
struct quern_mult_impl
{
    using unitary_data = quern_unitary_mat_mult<Val>;
    using Mat           = raw::Matrix<VM, Struct>;

    // form ud * mat
    static void eval_right(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        return quern_mult_struct<Val,VM,Struct>::eval_right(ret,mat,ud,t_unit);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        return quern_mult_struct<Val,VM,Struct>::eval_left(ret,mat,ud,t_unit);
    };
};

template<class Val, class Struct>
struct quern_mult_impl<Val, Integer, Struct>
{
    using unitary_data = quern_unitary_mat_mult<Val>;
    using VM            = Integer;
    using Mat           = raw::Matrix<VM, Struct>;
    using real_mat      = raw::Matrix<Real, Struct>;

    static void eval_right(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        real_mat rm = raw::converter<real_mat,Mat>::eval(mat);
        return quern_mult_impl<Val,Real,Struct>::eval_right(ret, rm, ud, t_unit);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const unitary_data& ud,
                           trans_type t_unit)
    {
        real_mat rm = raw::converter<real_mat,Mat>::eval(mat);
        return quern_mult_impl<Val,Real,Struct>::eval_left(ret, rm, ud, t_unit);
    };
};

//---------------------------------------------------------------
//              quern_unitary_mat_mult_right
//---------------------------------------------------------------
template<class Val>
struct quern_unitary_mat_mult_right : public extract_type_switch<void, quern_unitary_mat_mult_right<Val>,true>
{
    using unitary_data = quern_unitary_mat_mult<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, const unitary_data& ud, 
                     trans_type t_unit)
    {
        using VT = typename T::value_type;
        using ST = typename T::struct_type;
        return quern_mult_impl<Val,VT,ST>::eval_right(ret,mat,ud,t_unit);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, const unitary_data& ud, 
                            trans_type t_unit)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, ud, t_unit);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const unitary_data&, 
                     trans_type)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    }

    static void eval_scalar(const Matrix&, const Object&, Matrix&, const unitary_data&, 
                     trans_type)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    }
};

//---------------------------------------------------------------
//              quern_unitary_mat_mult_left
//---------------------------------------------------------------
template<class Val>
struct quern_unitary_mat_mult_left : public extract_type_switch<void, quern_unitary_mat_mult_left<Val>,true>
{
    using unitary_data = quern_unitary_mat_mult<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, const unitary_data& ud, 
                     trans_type t_unit)
    {
        using VT = typename T::value_type;
        using ST = typename T::struct_type;
        return quern_mult_impl<Val,VT,ST>::eval_left(ret,mat,ud,t_unit);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, const unitary_data& ud, 
                            trans_type t_unit)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, ud, t_unit);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const unitary_data&, 
                     trans_type)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    }

    static void eval_scalar(const Matrix&, const Object&, Matrix&, const unitary_data&, 
                     trans_type)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    }
};

//---------------------------------------------------------------
//              quern_unitary_mat_mult
//---------------------------------------------------------------
template<class Val>
quern_unitary_mat_mult<Val>::quern_unitary_mat_mult(Integer rows, Integer cols, const int* row_ind,
                                const int* col_ind, const Val* values)
    : m_rows(rows), m_cols(cols), row_start(row_ind), column_index(col_ind), value(values)
{};

template<class Val>
void quern_unitary_mat_mult<Val>::mult_right(Matrix& ret, const Matrix& X, 
                                             trans_type t_unit) const
{
    return quern_unitary_mat_mult_right<Val>::make<const Matrix&>(X,ret, *this, t_unit);
};
template<class Val>
void quern_unitary_mat_mult<Val>::mult_left(Matrix& ret, const Matrix& X, 
                                             trans_type t_unit) const
{
    return quern_unitary_mat_mult_left<Val>::make<const Matrix&>(X,ret, *this, t_unit);
};

template struct quern_unitary_mat_mult<Real>;
template struct quern_unitary_mat_mult<Float>;
template struct quern_unitary_mat_mult<Complex>;
template struct quern_unitary_mat_mult<Float_complex>;

}};