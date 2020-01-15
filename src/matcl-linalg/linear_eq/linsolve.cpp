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

#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "linsolve_utils.inl"
#include "linsolve_general.inl"
#include "linsolve_triang.inl"
#include "linsolve_qtriang.inl"
#include "linsolve_objects_decomp.h"
#include "linsolve_objects_ext.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/decompositions/chol.h"
#include "matcl-linalg/decompositions/balancing.h"

namespace matcl { namespace details
{

//--------------------------------------------------------------------------------
//                  VALUE TYPES
//--------------------------------------------------------------------------------
template<class S1, class S2, class V1, class V2, bool iso>
struct linsolve_general_impl
{
    using M1    = raw::Matrix<V1,S1>;
    using M2    = raw::Matrix<V2,S2>;

    static_assert(std::is_same<V1,Integer>::value == false
                  && std::is_same<V2,Integer>::value == false
                  && std::is_same<V1,Object>::value == false
                  && std::is_same<V2,Object>::value == false
                  , "invalid value types");

    static void eval(matcl::Matrix& ret, const M1& A, const permvec& p, 
                     const permvec& q, const M2& B, trans_type trans,
                     const matcl::options& opts)
    {
        return linsolve_general_struct<S1,S2,V1,V2>::eval(ret,A,p,q, B,trans,opts);
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, const permvec& p, 
                         const permvec& q, const M2& B, trans_type trans,
                         const matcl::options& opts)
    {
        return linsolve_general_struct<S1,S2,V1,V2>::eval_rev(ret,A,p,q, B,trans,opts);
    };
};

template<class S1, class S2, class V1, class V2>
struct linsolve_general_impl<S1,S2,V1,V2,true>
{
    using M1    = raw::Matrix<V1,S1>;
    using M2    = raw::Matrix<V2,S2>;

    static void eval(matcl::Matrix&, const M1&, const permvec&, const permvec&, 
                     const M2&, trans_type, const matcl::options&)
    {
        throw error::object_value_type_not_allowed("linsolve");
    };
    static void eval_rev(matcl::Matrix&, const M1&, const permvec&, const permvec&, 
                         const M2&, trans_type, const matcl::options&)
    {
        throw error::object_value_type_not_allowed("linsolve");
    };
};

template<class S1, class S2>
struct linsolve_general_impl<S1,S2,Integer,Integer,false>
{
    using M1    = raw::Matrix<Integer,S1>;
    using M2    = raw::Matrix<Integer,S2>;

    static void eval(matcl::Matrix& ret, const M1& A, const permvec& p, const permvec& q, 
                     const M2& B, trans_type trans, const matcl::options& opts)
    {
        using MT1   = raw::Matrix<Real,S1>;
        using MT2   = raw::Matrix<Real,S2>;

        MT1 Ac  = raw::converter<MT1,M1>::eval(A);
        MT2 Bc  = raw::converter<MT2,M2>::eval(B);

        return linsolve_general_impl<S1,S2,Real,Real,false>::eval(ret, Ac,p,q, Bc, trans, opts);
    };
    static void eval_rev(matcl::Matrix& ret, const M1& A, const permvec& p, const permvec& q, 
                         const M2& B, trans_type trans, const matcl::options& opts)
    {
        using MT1   = raw::Matrix<Real,S1>;
        using MT2   = raw::Matrix<Real,S2>;

        MT1 Ac  = raw::converter<MT1,M1>::eval(A);
        MT2 Bc  = raw::converter<MT2,M2>::eval(B);

        return linsolve_general_impl<S1,S2,Real,Real,false>::eval_rev(ret, Ac,p,q, Bc, trans, opts);
    };
};

template<class S, class V>
struct linsolve_general_lu_impl
{
    using M1    = raw::Matrix<V,S>;

    static void eval(linsolve_obj& ret, const M1& A, const matcl::options& opts)
    {
        return linsolve_general_lu_struct<S,V>::eval(ret,A,opts);
    };
};

template<class S>
struct linsolve_general_lu_impl<S, Integer>
{
    using M1    = raw::Matrix<Integer,S>;

    static void eval(linsolve_obj& ret, const M1& A, const matcl::options& opts)
    {
        using MT1   = raw::Matrix<Real,S>;

        MT1 Ac  = raw::converter<MT1,M1>::eval(A);

        return linsolve_general_lu_impl<S,Real>::eval(ret, Ac, opts);
    };
};

template<class S1, class S2, class V1, class V2, bool iso>
struct linsolve_triang_impl
{
    using M1    = raw::Matrix<V1,S1>;
    using M2    = raw::Matrix<V2,S2>;

    static_assert(std::is_same<V1,Integer>::value == false
                  && std::is_same<V2,Integer>::value == false
                  && std::is_same<V1,Object>::value == false
                  && std::is_same<V2,Object>::value == false
                  , "invalid value types");

    static void eval(matcl::Matrix& ret, const M1& A, const permvec& p, 
                     const permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const matcl::options& opts)
    {
        return linsolve_triang_str<S1,S2,V1,V2>::eval(ret, A,p,q, B,is_lt, trans, opts);
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, const permvec& p, 
                         const permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const matcl::options& opts)
    {
        return linsolve_triang_str<S1,S2,V1,V2>::eval_rev(ret, A,p,q, B,is_lt, trans, opts);
    };
};

template<class V, class S>
struct make_nonsing_str
{};

//-----------------------------------------------------------------
//                  make_nonsingular
//-----------------------------------------------------------------
template<class V>
struct make_nonsing_str<V,struct_dense>
{
    using ret_type  = tuple<Matrix,bool,bool>;
    using Mat       = raw::Matrix<V,struct_dense>;
    using VR        = typename md::real_type<V>::type;

    static void eval(ret_type& ret, const Mat& A, Real tol)
    {        
        const V* ptr    = A.ptr();
        Integer M       = A.rows();
        Integer N       = A.cols();
        Integer K       = std::min(M, N);
        Integer ld      = A.ld();

        if (tol == 0.0)
        {
            for (Integer i = 0; i < K; ++i)
            {
                if (mrd::is_zero(ptr[0]) == true)
                {
                    ret = ret_type(Matrix(A,false), true, false);
                    return;
                }

                ptr += ld + 1;
            }

            ret = ret_type(Matrix(A,false),false, false);
            return;
        };

        // calculate small pivot tolerance
        if (tol > 0.0)
        {
            Real max_val    = norm_vec_all(Matrix(A,false), basic_vector_norm::norm_inf);
            if (max_val == 0.0)
            {
                ret = ret_type(Matrix(A,false), true, false);
                return;
            }
            tol         = pow(constants::eps<VR>(), tol) * max_val;
        }
        else
        {
            tol         = -tol;
        }

        bool need_modif = false;
        VR tol_t        = VR(tol);
        Integer i;

        for (i = 0; i < K; ++i)
        {
            VR val  = abs(ptr[0]);

            if (val < tol_t)
            {
                need_modif = true;
                break;
            }

            ptr += ld + 1;
        }

        if (need_modif == false)
        {
            ret = ret_type(Matrix(A,false), false, false);
            return;
        };

        Mat A2  = A.make_unique();
        ld      = A2.ld();
        V* ptr2 = A2.ptr() + i * ld;        

        for (; i < K; ++i)
        {
            VR val  = abs(ptr2[0]);

            if (val < tol_t)
                ptr2[0]  = tol_t * (ptr2[0] == V(0) ? VR(1.0) : sign(ptr2[0]));

            ptr2 += ld + 1;
        }

        ret = ret_type(Matrix(A2,true), false, true);
    };
};

template<class V>
struct make_nonsing_str<V,struct_banded>
{
    using ret_type  = tuple<Matrix,bool,bool>;
    using Mat       = raw::Matrix<V,struct_banded>;
    using VR        = typename md::real_type<V>::type;

    static void eval(ret_type& ret, const Mat& A, Real tol)
    {    
        if (A.rows() == 0)
        {
            ret = ret_type(Matrix(A,false), false, false);
            return;
        };

        if (tol == 0.0)
        {
            if (A.has_diag(0) == false)
            {
                ret = ret_type(Matrix(A,false), true, false);
                return;
            }

            Integer s       = A.diag_length(0);
            const V* ptr    = A.rep_ptr() + A.first_elem_diag(0);
            Integer ld      = A.ld();

            for (Integer i = 0; i < s; ++i)
            {
                if (mrd::is_zero(ptr[0]) == true)
                {
                    ret = ret_type(Matrix(A,false), true, false);
                    return;
                }

                ptr += ld;
            }

            ret = ret_type(Matrix(A,false),false, false);
            return;
        };

        if (A.has_diag(0) == false)
        {
            Matrix B(A,false);
            B.diag(0)   = V(0.0);
            ret         = make_nonsingular(std::move(B), tol);
            return;
        };

        if (tol > 0.0)
        {
            Real max_val    = norm_vec_all(Matrix(A,false), basic_vector_norm::norm_inf);
            if (max_val == 0.0)
            {
                ret = ret_type(Matrix(A,false), true, false);
                return;
            }
            tol         = pow(constants::eps<VR>(), tol) * max_val;
        }
        else
        {
            tol         = -tol;
        };

        VR tol_t        = (VR)tol;
        bool need_modif = false;

        Integer i;
        {
            Integer s       = A.diag_length(0);
            const V* ptr    = A.rep_ptr() + A.first_elem_diag(0);
            Integer ld      = A.ld();

            for (i = 0; i < s; ++i)
            {
                VR val      = abs(ptr[0]);
                if (val < tol_t)
                {
                    need_modif  = true;
                    break;
                }

                ptr += ld;
            }
        };

        if (need_modif == false)
        {
            ret = ret_type(Matrix(A,false), false, false);
            return;
        };

        Mat A2          = A.make_unique();

        {
            Integer s       = A2.diag_length(0);
            V* ptr          = A2.rep_ptr() + A2.first_elem_diag(0);
            Integer ld      = A2.ld();

            for (; i < s; ++i)
            {
                VR val  = abs(ptr[0]);

                if (val < tol_t)
                    ptr[0]  = tol_t * (ptr[0] == V(0) ? VR(1.0) : sign(ptr[0]));


                ptr += ld;
            }
        };

        ret = ret_type(Matrix(A2,true), false, true);
    };
};

template<class V>
struct make_nonsing_str<V,struct_sparse>
{
    using ret_type  = tuple<Matrix,bool,bool>;
    using Mat       = raw::Matrix<V,struct_sparse>;
    using VR        = typename md::real_type<V>::type;

    static void eval(ret_type& ret, const Mat& A, Real tol)
    {     
        Integer N                   = A.rows();

        if (N == 0)
        {
            ret = ret_type(Matrix(A,false), false, false);
            return;
        };

        if (A.nnz() == 0)
        {
            value_code vc   = matrix_traits::value_code<VR>::value;

            if (tol >= 0.0)
            {
                ret = ret_type(Matrix(A,false), true, false);
                return;
            }
            else
            {
                tol = -tol;
                ret = ret_type(VR(tol) * speye(N,N,vc), false, true);
                return;
            }
        };

        using ccs   = mrd::sparse_ccs<V>;

        if (tol == 0.0)
        {
            const ccs& rep          = A.rep();
            const V* ptr_x          = rep.ptr_x();

            for (Integer i = 0; i < N; ++i)
            {
                Integer k;
                bool has    = rep.has_element(i,i,k);

                if (has == false || mrd::is_zero(ptr_x[k]))
                {
                    ret = ret_type(Matrix(A,false), true, false);
                    return;
                }
            }

            ret = ret_type(Matrix(A,false),false, false);
            return;
        };

        // calculate small pivot tolerance
        if (tol > 0.0)
        {
            Real max_val    = norm_vec_all(Matrix(A,false), basic_vector_norm::norm_inf);
            if (max_val == 0.0)
            {
                ret = ret_type(Matrix(A,false), true, false);
                return;
            }

            tol         = pow(constants::eps<VR>(), tol) * max_val;
        }
        else
        {
            tol         = -tol;
        };                    

        VR tol_t        = (VR)tol;
        bool need_modif = false;
        Integer i;

        {
            const ccs& rep  = A.rep();
            const V* ptr_x  = rep.ptr_x();

            for (i = 0; i < N; ++i)
            {
                Integer k;
                bool has    = rep.has_element(i,i,k);

                if (has == false || abs(ptr_x[k]) < tol_t)
                {
                    need_modif  = true;
                    break;
                }
            }
        };

        if (need_modif == false)
        {
            ret = ret_type(Matrix(A,false), false, false);
            return;
        };

        Matrix B            = matcl::convert(Matrix(A, false), Mat::matrix_code);
        B.diag(0).add_sparse();
        Mat& A2             = B.get_impl_unique<Mat>();

        {
            ccs& rep        = A2.rep();
            V* ptr_x        = rep.ptr_x();

            for (; i < N; ++i)
            {
                Integer k;
                bool has    = rep.has_element(i,i,k);

                matcl_assert(has == true, "invalid sparse structure");

                VR val          = abs(ptr_x[k]);

                if (val < tol_t)
                    ptr_x[k]    = tol_t * (ptr_x[k] == V(0) ? VR(1.0) : sign(ptr_x[k]));
            }
        };

        ret = ret_type(Matrix(A2,true), false, true);
    };
};

template<class V, class S>
struct make_nonsing_impl
{
    using ret_type  = tuple<Matrix,bool,bool>;
    using Mat       = raw::Matrix<V,S>;

    static void eval(ret_type& ret, const Mat& A, Real tol)
    {
        return make_nonsing_str<V,S>::eval(ret, A, tol);
    };
};
template<class S>
struct make_nonsing_impl<Integer,S>
{
    using ret_type  = tuple<Matrix,bool,bool>;
    using Mat       = raw::Matrix<Integer,S>;

    static void eval(ret_type& ret, const Mat& A, Real tol)
    {
        using Mat_R = raw::Matrix<Real,S>;
        Mat_R Ar    = raw::converter<Mat_R, Mat>::eval(A);
        return make_nonsing_impl<Real,S>::eval(ret, Ar, tol);
    };
};

//-----------------------------------------------------------------
//              linsolve_triang_impl
//-----------------------------------------------------------------
template<class S1, class S2, class V1, class V2>
struct linsolve_triang_impl<S1,S2,V1,V2,true>
{
    using M1    = raw::Matrix<V1,S1>;
    using M2    = raw::Matrix<V2,S2>;
    
    static void eval(matcl::Matrix&, const M1& , const permvec&, const permvec&, 
                     const M2& , bool, trans_type, const matcl::options&)
    {
        throw error::object_value_type_not_allowed("linsolve");
    };

    static void eval_rev(matcl::Matrix&, const M1&, const permvec&, const permvec&, 
                         const M2& , bool, trans_type, const matcl::options&)
    {
        throw error::object_value_type_not_allowed("linsolve");
    };
};

template<class S1, class S2>
struct linsolve_triang_impl<S1,S2,Integer,Integer,false>
{
    using M1    = raw::Matrix<Integer,S1>;
    using M2    = raw::Matrix<Integer,S2>;

    static void eval(matcl::Matrix& ret, const M1& A, const permvec& p, const permvec& q, 
                     const M2& B, bool is_lt, trans_type trans, const matcl::options& opts)
    {
        using MT1   = raw::Matrix<Real,S1>;
        using MT2   = raw::Matrix<Real,S2>;

        MT1 Ac  = raw::converter<MT1,M1>::eval(A);
        MT2 Bc  = raw::converter<MT2,M2>::eval(B);

        return linsolve_triang_impl<S1,S2,Real,Real,false>::eval(ret, Ac,p,q, Bc,is_lt,trans, opts);
    };
    static void eval_rev(matcl::Matrix& ret, const M1& A, const permvec& p, const permvec& q, 
                         const M2& B, bool is_lt, trans_type trans, const matcl::options& opts)
    {
        using MT1   = raw::Matrix<Real,S1>;
        using MT2   = raw::Matrix<Real,S2>;

        MT1 Ac  = raw::converter<MT1,M1>::eval(A);
        MT2 Bc  = raw::converter<MT2,M2>::eval(B);

        return linsolve_triang_impl<S1,S2,Real,Real,false>::eval_rev(ret, Ac,p,q, Bc,is_lt,trans, opts);
    };
};

template<class S1, class S2, class V1, class V2, bool iso>
struct linsolve_qtriang_impl
{
    using M1    = raw::Matrix<V1,S1>;
    using M2    = raw::Matrix<V2,S2>;

    static_assert(std::is_same<V1,Integer>::value == false
                  && std::is_same<V2,Integer>::value == false
                  && std::is_same<V1,Object>::value == false
                  && std::is_same<V2,Object>::value == false
                  , "invalid value types");

    static void eval(matcl::Matrix& ret, const M1& A, const permvec& p, 
                     const permvec& q, const M2& B, trans_type trans,
                     const matcl::options& opts)
    {
        return linsolve_qtriang_str<S1,S2,V1,V2>::eval(ret, A,p,q, B, trans, opts);
    };
};

template<class S1, class S2, class V1, class V2>
struct linsolve_qtriang_impl<S1, S2, V1, V2, true>
{
    using M1    = raw::Matrix<V1,S1>;
    using M2    = raw::Matrix<V2,S2>;

    static void eval(matcl::Matrix&, const M1&, const permvec&, 
                     const permvec&, const M2&, trans_type,
                     const matcl::options&)
    {
        throw error::object_value_type_not_allowed("linsolve");
    };
};

template<class S1, class S2>
struct linsolve_qtriang_impl<S1,S2,Integer,Integer,false>
{
    using M1    = raw::Matrix<Integer,S1>;
    using M2    = raw::Matrix<Integer,S2>;

    static void eval(matcl::Matrix& ret, const M1& A, const permvec& p, const permvec& q, 
                     const M2& B, trans_type trans, const matcl::options& opts)
    {
        using MT1   = raw::Matrix<Real,S1>;
        using MT2   = raw::Matrix<Real,S2>;

        MT1 Ac  = raw::converter<MT1,M1>::eval(A);
        MT2 Bc  = raw::converter<MT2,M2>::eval(B);

        return linsolve_qtriang_impl<S1,S2,Real,Real,false>::eval(ret, Ac,p,q, Bc, trans, opts);
    };
};

//--------------------------------------------------------------------------------
//                  SWITCH FUNCTORS
//--------------------------------------------------------------------------------

struct linsolve_triang_functor : public extract_type2_switch<void,linsolve_triang_functor,
                                        matcl::raw::val_type_corrector_int, ver_nonstatic>
{
    permvec         m_p;
    permvec         m_q;
    trans_type      m_trans;
    bool            is_lt;
    const options&  m_opts;

    linsolve_triang_functor(const matcl::permvec& p, const matcl::permvec& q, 
                            bool is_lt, trans_type trans, const options& opts) 
        : is_lt(is_lt), m_trans(trans), m_p(p), m_q(q), m_opts(opts)
    {};

    template<class MT1, class MT2>
    void eval_mat_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using S1    = typename MT1::struct_type;
        using S2    = typename MT2::struct_type;
        using V1    = typename MT1::value_type;
        using V2    = typename MT2::value_type;

        static const bool iso = std::is_same<V1,Object>::value ||
                                std::is_same<V2,Object>::value;

        return linsolve_triang_impl<S1,S2,V1,V2,iso>
            ::eval(ret, A,m_p,m_q, B, is_lt, m_trans, m_opts);
    };

    template<class MT1, class MT2>
    void eval_mat_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM    = raw::Matrix<MT2,struct_dense>;
        DM Bc(ti::get_ti(B),B,1,1);

        return eval_mat_mat<MT1,DM>(A,Bc, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM    = raw::Matrix<MT1,struct_dense>;
        DM Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<DM,MT2>(Ac,B, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        if (A == MT1(0.))
            throw error::error_singular();

        typename correct_type<MT1>::type Ac = correct_type<MT1>::eval(A);
        typename correct_type<MT2>::type Bc = correct_type<MT2>::eval(B);
        
        if (this->m_trans == trans_type::conj_trans)
            ret = Bc/matcl::conj(Ac);
        else
            ret = Bc/Ac;
    };
};

struct linsolve_triang_functor_rev : public extract_type2_switch<void,linsolve_triang_functor_rev,
                                        matcl::raw::val_type_corrector_int, ver_nonstatic>
{
    permvec         m_p;
    permvec         m_q;
    trans_type      m_trans;
    bool            is_lt;
    const options&  m_opts;

    linsolve_triang_functor_rev(const matcl::permvec& p, const matcl::permvec& q, 
                                bool is_lt, trans_type trans, const options& opts) 
        : is_lt(is_lt), m_trans(trans), m_p(p), m_q(q), m_opts(opts)
    {};

    template<class MT1, class MT2>
    void eval_mat_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using S1    = typename MT1::struct_type;
        using S2    = typename MT2::struct_type;
        using V1    = typename MT1::value_type;
        using V2    = typename MT2::value_type;

        static const bool iso = std::is_same<V1,Object>::value ||
                                std::is_same<V2,Object>::value;

        return linsolve_triang_impl<S1,S2,V1,V2,iso>
                ::eval_rev(ret, A,m_p,m_q, B,is_lt, m_trans, m_opts);
    };

    template<class MT1, class MT2>
    void eval_mat_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT2,struct_dense>;
        DM Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<MT1,DM>(A,Bc, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT1,struct_dense>;
        DM Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<DM,MT2>(Ac,B, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        typename correct_type<MT1>::type Ac = correct_type<MT1>::eval(A);
        typename correct_type<MT2>::type Bc = correct_type<MT2>::eval(B);

        if (Ac == 0.)
            throw error::error_singular();
        
        if (this->m_trans == trans_type::conj_trans)
            ret = Bc/matcl::conj(Ac);
        else
            ret = Bc/Ac;
    };
};

struct linsolve_qtriang_functor : public extract_type2_switch<void,linsolve_qtriang_functor,
                                    matcl::raw::val_type_corrector_int,ver_nonstatic>
{
    permvec         m_p;
    permvec         m_q;
    trans_type      m_trans;
    const options&  m_opts;

    linsolve_qtriang_functor(const matcl::permvec& p, const matcl::permvec& q, 
                            trans_type trans, const options& opts) 
        : m_trans(trans), m_p(p), m_q(q), m_opts(opts)
    {};

    template<class MT1, class MT2>
    void eval_mat_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using S1    = typename MT1::struct_type;
        using S2    = typename MT2::struct_type;
        using V1    = typename MT1::value_type;
        using V2    = typename MT2::value_type;

        static const bool iso = std::is_same<V1,Object>::value ||
                                std::is_same<V2,Object>::value;

        return linsolve_qtriang_impl<S1,S2,V1,V2,iso>
            ::eval(ret, A, m_p,m_q, B, m_trans, m_opts);
    };

    template<class MT1, class MT2>
    void eval_mat_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT2,struct_dense>;
        DM Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<MT1,DM>(A,Bc, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT1,struct_dense>;
        DM Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<DM,MT2>(Ac,B, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        typename correct_type<MT1>::type Ac = correct_type<MT1>::eval(A);
        typename correct_type<MT2>::type Bc = correct_type<MT2>::eval(B);

        if (Ac == 0.)
            throw error::error_singular();
        
        if (this->m_trans == trans_type::conj_trans)
            ret = Bc/matcl::conj(Ac);
        else
            ret = Bc/Ac;
    };
};

struct linsolve_general_functor : public extract_type2_switch<void,linsolve_general_functor,
                                    matcl::raw::val_type_corrector_int, ver_nonstatic>
{
    permvec         m_p;
    permvec         m_q;
    trans_type      m_trans;
    const options&  m_opts;
    
    linsolve_general_functor(const permvec& p, const permvec& q, trans_type trans,
                             const options& opts)
        :m_trans(trans), m_p(p), m_q(q), m_opts(opts)
    {};

    template<class MT1, class MT2>
    void eval_mat_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using S1    = typename MT1::struct_type;
        using S2    = typename MT2::struct_type;
        using V1    = typename MT1::value_type;
        using V2    = typename MT2::value_type;

        static const bool iso = std::is_same<V1,Object>::value ||
                                std::is_same<V2,Object>::value;

        return linsolve_general_impl<S1,S2,V1,V2,iso>
            ::eval(ret, A,m_p,m_q, B, m_trans, m_opts);
    };

    template<class MT1, class MT2>
    void eval_mat_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT2,struct_dense>;
        DM Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<MT1,DM>(A,Bc, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT1,struct_dense>;
        DM Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<DM,MT2>(Ac,B, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        typename correct_type<MT1>::type Ac = correct_type<MT1>::eval(A);
        typename correct_type<MT2>::type Bc = correct_type<MT2>::eval(B);

        if (Ac == 0.)
            throw error::error_singular();
        
        if (this->m_trans == trans_type::conj_trans)
            ret =  Bc / matcl::conj(Ac);            
        else
            ret = Bc/Ac;
    };
};

struct linsolve_general_lu_functor : public extract_type_switch<void,linsolve_general_lu_functor,true>
{
    using linsolve_ptr  = linsolve_obj::linsolve_data_ptr;

    template<class T>
    static void eval(const Matrix& h, const T& mat, linsolve_obj& ret, const matcl::options& opts)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        if (mat.rows() == 1)
        {
            if (is_nan(mat(1,1)) == true)
                ret = linsolve_nan(1, matrix_traits::value_code<V>::value);
            else
                ret = linsolve_obj(linsolve_ptr(new details::linsolve_obj_scalar(h, opts)));                

            return;
        };

        return linsolve_general_lu_impl<S,V>::eval(ret, mat, opts);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix&, const T& A, linsolve_obj& ret, const matcl::options& opts)
    {
        if (is_nan(A) == true)
            ret = linsolve_nan(1, matrix_traits::value_code<T>::value);
        else
            ret = linsolve_obj(linsolve_ptr(new details::linsolve_obj_scalar(A, opts)));
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, linsolve_obj&, const matcl::options&)
    {
        throw error::object_value_type_not_allowed("linsolve_lu");
    };
    static void eval_scalar(const Matrix&, const Object&, linsolve_obj&, const matcl::options&)
    {
        throw error::object_value_type_not_allowed("linsolve_lu");
    };
};

struct vis_make_nonsingular : public extract_type_switch<void,vis_make_nonsingular,true>
{
    using ret_type = tuple<Matrix,bool,bool>;

    template<class T>
    static void eval(const Matrix&, const T& mat, ret_type& ret, Real tol)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return make_nonsing_impl<V, S>::eval(ret, mat, tol);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix&, const T& A, ret_type& ret, Real tol)
    {
        if (mrd::is_zero(A) == true && tol >= 0.0)
            ret = ret_type(A, true, false);

        using VR    = typename md::real_type<T>::type;
        tol         = -tol;

        if (abs(A) < tol)
            ret = ret_type(VR(tol) * sign(A), false, true);
        else
            ret = ret_type(A, false, false);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, ret_type&, Real)
    {
        throw error::object_value_type_not_allowed("make_nonsingular");
    };
    static void eval_scalar(const Matrix&, const Object&, ret_type&, Real)
    {
        throw error::object_value_type_not_allowed("make_nonsingular");
    };
};
    
struct linsolve_general_functor_rev : public extract_type2_switch<void,linsolve_general_functor_rev,
                                        matcl::raw::val_type_corrector_int, ver_nonstatic>
{
    permvec         m_p;
    permvec         m_q;
    trans_type      m_trans;
    const options&  m_opts;
    
    linsolve_general_functor_rev(const matcl::permvec& p, const matcl::permvec& q, 
                                 trans_type trans, const options& opts)
        :m_trans(trans), m_p(p), m_q(q), m_opts(opts)
    {};

    template<class MT1, class MT2>
    void eval_mat_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using S1    = typename MT1::struct_type;
        using S2    = typename MT2::struct_type;
        using V1    = typename MT1::value_type;
        using V2    = typename MT2::value_type;

        static const bool iso = std::is_same<V1,Object>::value ||
                                std::is_same<V2,Object>::value;

        return linsolve_general_impl<S1,S2,V1,V2,iso>::eval_rev(ret, A,m_p,m_q, B, m_trans, m_opts);
    };

    template<class MT1, class MT2>
    void eval_mat_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT2,struct_dense>;
        DM Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<MT1,DM>(A,Bc, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_mat(const MT1& A, const MT2& B, Matrix& ret)
    {
        using DM = raw::Matrix<MT1,struct_dense>;
        DM Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<DM,MT2>(Ac,B, ret);
    };

    template<class MT1, class MT2>
    void eval_scal_scal(const MT1& A, const MT2& B, Matrix& ret)
    {
        typename correct_type<MT1>::type Ac = correct_type<MT1>::eval(A);
        typename correct_type<MT2>::type Bc = correct_type<MT2>::eval(B);

        if (Ac == 0.)
            throw error::error_singular();
        
        if (this->m_trans == trans_type::conj_trans)
            ret = Bc / matcl::conj(Ac);
        else
            ret = Bc/Ac;
    };
};

static void linsolve_id(Matrix& ret, const matcl::permvec& p, const matcl::permvec& q, 
                   const matcl::Matrix& b, trans_type trans, const matcl::options& opts, value_code vc)
{
    (void)opts;

    if (p.is_id() == false && q.is_id() == false)
    {
        if (trans == trans_type::no_trans)
        {
            Matrix out = b(p(q.invperm()).to_matrix(),colon());
            ret = details::convert_value(out, vc);
            return;
        }
        else
        {
            Matrix out = b(q(p.invperm()).to_matrix(),colon());
            ret = details::convert_value(out, vc);
            return;
        };
    }
    else if (p.is_id() == false)
    {
        if (trans == trans_type::no_trans)
        {
            Matrix out = b(p.to_matrix(),colon());
            ret = details::convert_value(out, vc);
            return;
        }
        else
        {
            Matrix out = b(p.invperm().to_matrix(),colon());
            ret = details::convert_value(out, vc);
            return;
        };
    }
    else if (q.is_id() == false)
    {
        if (trans == trans_type::no_trans)
        {
            Matrix out = b(q.invperm().to_matrix(),colon());
            ret = details::convert_value(out, vc);
            return;
        }
        else
        {
            Matrix out = b(q.to_matrix(),colon());
            ret = details::convert_value(out, vc);
            return;
        };
    }
    else
    {
        ret = details::convert_value(b, vc);
        return;
    };
};

static void linsolve_id_rev(Matrix& ret, const matcl::permvec& p, const matcl::permvec& q, 
                       const matcl::Matrix& b, trans_type trans, const matcl::options& opts, value_code vc)
{
    (void)opts;

    if (p.is_id() == false && q.is_id() == false)
    {
        if (trans == trans_type::no_trans)
        {
            Matrix out = b(colon(), q(invperm(p)).to_matrix());
            ret = details::convert_value(out, vc);
            return;
        }
        else
        {
            Matrix out = b(colon(),p(invperm(q)).to_matrix());
            ret = details::convert_value(out, vc);
            return;
        };
    }
    else if (p.is_id() == false)
    {
        if (trans == trans_type::no_trans)
        {
            Matrix out = b(colon(), invperm(p).to_matrix());
            ret = details::convert_value(out, vc);
            return;
        }
        else
        {
            Matrix out = b(colon(),p.to_matrix());
            ret = details::convert_value(out, vc);
            return;
        };
    }
    else if (q.is_id() == false)
    {
        if (trans == trans_type::no_trans)
        {
            Matrix out = b(colon(), q.to_matrix());
            ret = details::convert_value(out, vc);
            return;
        }
        else
        {
            Matrix out = b(colon(),invperm(q).to_matrix());
            ret = details::convert_value(out, vc);
            return;
        };
    }
    else
    {
        ret = details::convert_value(b, vc);
        return;
    };
};

static void linsolve_diag(Matrix& ret, const matcl::Matrix& A, matcl::permvec p, 
                          matcl::permvec q, const matcl::Matrix& b, trans_type trans,
                          const matcl::options& opts)
{
    (void)opts;

    matcl::Matrix Ad        = get_diag(A);
    bool is_nsing           = all_vec(Ad);

    if (is_nsing == false)
        throw error::error_singular();

    Matrix out;
    Matrix b2 = b;

    if (trans != trans_type::no_trans)
        std::swap(p,q);

    if (p.is_id() == false)
        b2  = b2(p.to_matrix(),colon());

    value_code vc_A = Ad.get_value_code();
    value_code vc_B = b2.get_value_code();

    value_code vc   = matrix_traits::unify_value_types(vc_A, vc_B);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
    value_code vcr  = matrix_traits::real_value_type(vc);

    Matrix one      = ones(1, 1, vcr);

    switch(trans)
    {
        case trans_type::no_trans:
        case trans_type::trans:
        {
            out = bdiag(div(one, Ad)) * b2;
            break;
        }
        case trans_type::conj_trans:
        {
            Ad = matcl::conj(Ad);
            out = bdiag(div(one, Ad)) * b2;
            break;
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };    

    if (q.is_id() == false)
        out = out(q.invperm(), colon());

    ret = out;
    return;
};

static void linsolve_diag_rev(Matrix& ret, const matcl::Matrix& A, permvec p, permvec q,
                         const matcl::Matrix& b, trans_type trans, const matcl::options& opts)
{
    (void)opts;

    matcl::Matrix Ad    = get_diag(A);
    bool is_nsing       = all_vec(Ad);

    if (is_nsing == false)
        throw error::error_singular();

    Matrix b2 = b;

    if (trans != trans_type::no_trans)
        std::swap(p,q);

    if (q.is_id() == false)
        b2 = b2(colon(),q.to_matrix());

    Matrix out;

    value_code vc_A = Ad.get_value_code();
    value_code vc_B = b2.get_value_code();

    value_code vc   = matrix_traits::unify_value_types(vc_A, vc_B);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
    value_code vcr  = matrix_traits::real_value_type(vc);

    Matrix one      = ones(1, 1, vcr);

    switch(trans)
    {
        case trans_type::no_trans:
        case trans_type::trans:
        {
            out = b2 * bdiag(div(one, Ad));
            break;
        }
        case trans_type::conj_trans:
        {
            Ad = matcl::conj(Ad);
            out = b2 * bdiag(div(one, Ad));
            break;
        }
        default:
        {
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
        }
    };    

    if (p.is_id() == false)
        out = out(colon(),invperm(p).to_matrix());

    ret = out;
    return;
};

static void linsolve_unitary(Matrix& ret, const matcl::Matrix& A, matcl::permvec p, 
                             matcl::permvec q, const matcl::Matrix& b, trans_type trans,
                             const matcl::options& opts)
{
    (void)opts;

    Matrix b2 = b;

    if (trans != trans_type::no_trans)
        std::swap(p,q);

    if (p.is_id() == false)
        b2  = b2(p.to_matrix(),colon());

    Matrix out;

    switch(trans)
    {
        case trans_type::no_trans:
            out = matcl::ctrans(A)*b2;
            break;
        case trans_type::trans:
            out = matcl::conj(A)*b2;
            break;
        case trans_type::conj_trans:
            out = A*b2;
            break;
        default:
        {
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
        }
    };     

    if (q.is_id() == false)
        out  = out(invperm(q).to_matrix(),colon());

    ret = out;
    return;
};

static void linsolve_unitary_rev(Matrix& ret, const matcl::Matrix& A, permvec p, 
                                        permvec q, const matcl::Matrix& b, trans_type trans,
                                        const matcl::options& opts)
{
    (void)opts;

    Matrix b2 = b;

    if (trans != trans_type::no_trans)
        std::swap(p,q);

    if (q.is_id() == false)
        b2 = b2(colon(),q.to_matrix());

    Matrix out;

    switch(trans)
    {
        case trans_type::no_trans:
            out = b2*matcl::ctrans(A);
            break;
        case trans_type::trans:
            out = b2*matcl::conj(A);
            break;
        case trans_type::conj_trans:
            out = b2*A;
            break;
        default:
        {
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
        }
    };     

    if (p.is_id() == false)
        out = out(colon(),invperm(p).to_matrix());

    ret = out;
    return;
};

static void linsolve_triang(Matrix& ret, const matcl::Matrix& A, const matcl::permvec& p, 
                            const matcl::permvec& q, const matcl::Matrix& b,bool is_lt, trans_type trans,
                            const matcl::options& opts)
{
    return details::linsolve_triang_functor(p,q,is_lt,trans, opts).make(A,b, ret);
};

static void linsolve_triang_rev(Matrix& ret, const matcl::Matrix& A, const permvec& p, 
                                const permvec& q, const matcl::Matrix& b,bool is_lt, trans_type trans,
                                const matcl::options& opts)
{
    return details::linsolve_triang_functor_rev(p,q,is_lt,trans, opts).make(A,b, ret);
};

static void linsolve_qtriang(Matrix& ret, const matcl::Matrix& A, const matcl::permvec& p, 
                             const matcl::permvec& q,  const matcl::Matrix& b,trans_type trans,
                             const matcl::options& opts)
{
    return details::linsolve_qtriang_functor(p,q,trans, opts).make(A,b, ret);
};

static void linsolve_general(Matrix& ret, const matcl::Matrix& A, const permvec& p, 
                             const permvec& q, const matcl::Matrix& b, trans_type trans,
                             const matcl::options& opts)
{
    return details::linsolve_general_functor(p,q,trans,opts).make(A,b,ret);
};

static void linsolve_lu_general(linsolve_obj& ret, const matcl::Matrix& A, const matcl::options& opts)
{
    return details::linsolve_general_lu_functor::make<const Matrix&>(A,ret,opts);
};

static void linsolve_general_rev(Matrix& ret, const matcl::Matrix& A, const permvec& p, 
                                 const permvec& q, const matcl::Matrix& b, trans_type trans,
                                 const matcl::options& opts)
{
    return details::linsolve_general_functor_rev(p,q,trans, opts).make(A,b, ret);
};

static void linsolve_impl(Matrix& ret, const matcl::Matrix& A, const matcl::permvec& p, 
                const matcl::permvec& q, const matcl::Matrix& b, trans_type trans,
                const matcl::options& opts)
{
    //A must be square
    error::check_lsolve(A.rows(), A.cols(), b.rows(), b.cols());

    //p,q must be permutation vectors
    if (p.length() != A.rows())
        throw error::invalid_permvec_length(p.length(), A.rows());

    if (q.length() != A.cols())
        throw error::invalid_permvec_length(q.length(), A.cols());

    Integer N = A.rows();

    value_code vc_A = A.get_value_code();
    value_code vc_B = b.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_A, vc_B);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
    value_code vcr  = matrix_traits::real_value_type(vc);

    if (vc == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve");

    if (b.cols() == 0 || N == 0)
    {
        ret = details::convert_value(Matrix(b), vc);
        return;
    }

    if (b.structural_nnz() == 0)
    {
        ret = details::convert_value(Matrix(b), vc);
        return;
    };

    if (N == 1)
    {
        if (A == 0)
            throw error::error_singular();

        Matrix one  = ones(1,1, vcr);

        if (trans == trans_type::conj_trans)
        {
            ret = b * div(one, matcl::conj(A));
            return;
        }
        else
        {
            ret = b * div(one, A);
            return;
        };
    }

    if (A.get_struct().is_id())
        return linsolve_id(ret,p,q,b, trans,opts, vc);

    Integer ldiags      = matcl::get_ld(A, 0);
    Integer udiags      = matcl::get_ud(A, 0);

    using pivot_type    = opt::linsolve::pivot_type;
    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;
    bool isp            = (piv == pivot_type::partial);

    if (ldiags == 0 && udiags == 0)
        return linsolve_diag(ret,A,p,q, b, trans,opts);

    if (isp && ldiags == 0)
        return linsolve_triang(ret,A,p,q, b,false, trans,opts);

    if (isp && udiags == 0)
        return linsolve_triang(ret,A,p,q, b,true, trans, opts);

    if (is_unitary(A.get_struct()))
        return linsolve_unitary(ret,A,p,q, b, trans,opts);

    //now quasi-triangular matrices
    if (isp && A.get_struct().is_qtriu())
    {
        bool is_real_A      = matrix_traits::is_float(vc_A) || vc_A == value_code::v_integer;
        bool is_real_B      = matrix_traits::is_float(vc_B) || vc_B == value_code::v_integer;

        //implemented only for dense matrices real matrices
        if (A.get_struct_code() == struct_code::struct_dense && is_real_A && is_real_B)
            return linsolve_qtriang(ret, A,p,q, b, trans, opts);
    };

    return linsolve_general(ret, A,p,q, b, trans, opts);
};

static void linsolve_lu_impl(linsolve_obj& ret, const matcl::Matrix& A, const matcl::options& opts)
{
    //A must be square
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    using linsolve_ptr  = linsolve_obj::linsolve_data_ptr;

    value_code vc_A = A.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_A, value_code::v_float);

    if (vc == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_lu");

    Integer ldiags      = matcl::get_ld(A, 0);
    Integer udiags      = matcl::get_ud(A, 0);

    if (ldiags == 0 && udiags == 0)
    {
        ret = linsolve_diag(A, opts);
        return;
    }

    if (is_unitary(A.get_struct()))
    {
        ret = linsolve_unitary(A);
        return;
    }

    using pivot_type    = opt::linsolve::pivot_type;
    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;
    bool isp            = (piv == pivot_type::partial);

    if (isp && udiags == 0)
    {
        ret = linsolve_triang(A, opts);
        return;
    }

    if (isp && ldiags == 0)
    {
        ret = linsolve_triang(A, opts);
        return;
    }

    //now quasi-triangular matrices
    if (isp && A.get_struct().is_qtriu())
    {
        ret = linsolve_uhess(A, opts);
        return;
    };
    if (isp && A.get_struct().is_qtril() && A.get_struct_code() == struct_code::struct_dense)
    {
        ret = linsolve_lhess(A, opts);
        return;
    };

    if (isp && ldiags == 1 && udiags == 1)
    {
        Matrix DL = A.diag(-1);
        Matrix D0 = A.diag(0);        
        Matrix DU = A.diag(1);

        ret = linsolve_tridiag(DL, D0, DU, opts);
        return;
    };

    bool balance    = false;
    Real tol_sig    = 0.0;

    if ((isp == true) && opts.get_option<bool>(opt::linsolve::do_balancing()))
    {
        balance     = true;
        tol_sig     = 0.0;
    }
    else if ((isp == false) && opts.get_option<bool>(opt::linsolve::do_balancing_rr()))
    {
        balance     = true;

        tol_sig     = opts.get_option<Real>(opt::linsolve::tol_sing());
        if (tol_sig < 0.0)
        {
            tol_sig = - tol_sig;
        }
        else if (tol_sig > 0.0)
        {
            tol_sig = norm(A, -2.0) * pow(constants::eps(A.get_value_code()), tol_sig);
        };
    };

    if (balance == true)
    {
        Matrix B, R, C;
        tie(B,R,C) = balance_gen2(A, true, tol_sig);

        linsolve_obj lo_B;
        linsolve_lu_general(lo_B, B, opts);

        ret = linsolve_balanced(A, R, lo_B, C);
        return;
    }
    else
    {
        return linsolve_lu_general(ret, A, opts);
    };
};

static void linsolve_rev_impl(Matrix& ret, const matcl::Matrix& A, const matcl::permvec& p, 
                              const matcl::permvec& q, const matcl::Matrix& b,
                              const matcl::options& opts)
{
    //transposed versions cannot be efficiently implemented for sparse A
    //matrix

    //A must be square
    error::check_lsolve_rev(A.rows(), A.cols(), b.rows(), b.cols());
    
    //p,q must be permutation vectors
    if (p.length() != A.rows())
        throw error::invalid_permvec_length(p.length(), A.rows());

    if (q.length() != A.cols())
        throw error::invalid_permvec_length(q.length(), A.cols());

    Integer N = A.rows();

    value_code vc_A = A.get_value_code();
    value_code vc_B = b.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_A, vc_B);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
    value_code vcr  = matrix_traits::real_value_type(vc);

    if (vc == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_rev");

    if (b.rows() == 0 || N == 0)
    {
        ret = details::convert_value(Matrix(b), vc);
        return;
    }

    if (b.structural_nnz() == 0)
    {
        ret = details::convert_value(Matrix(b), vc);
        return;
    }

    if (N == 1)
    {
        if (A == 0)
            throw error::error_singular();

        Matrix one  = ones(1,1, vcr);
        ret         = b * div(one, A);
        return;
    }

    if (A.get_struct().is_id())
        return linsolve_id_rev(ret, p,q,b, trans_type::no_trans, opts, vc);

    Integer ldiags      = matcl::get_ld(A, 0);
    Integer udiags      = matcl::get_ud(A, 0);

    using pivot_type    = opt::linsolve::pivot_type;
    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;
    bool isp            = (piv == pivot_type::partial);

    if (ldiags == 0 && udiags == 0)
        return linsolve_diag_rev(ret, A,p,q, b, trans_type::no_trans, opts);

    if (isp && ldiags == 0)
        return linsolve_triang_rev(ret, A,p,q, b,false, trans_type::no_trans, opts);

    if (isp && udiags == 0)
        return linsolve_triang_rev(ret, A,p,q, b,true, trans_type::no_trans, opts);

    if (is_unitary(A.get_struct()))
        return linsolve_unitary_rev(ret, A,p,q, b, trans_type::no_trans, opts);

    //quasi-triangular matrices case is not implemented

    return linsolve_general_rev(ret, A,p,q, b, trans_type::no_trans, opts);
};

static void linsolve_rev2_impl(Matrix& ret, const matcl::Matrix& A, const matcl::permvec& p, 
                            const matcl::permvec& q, const matcl::Matrix& b, trans_type trans,
                            const matcl::options& opts)
{
    //transposed versions cannot be efficiently implemented for sparse A
    //matrix

    //A must be square
    error::check_lsolve_rev(A.rows(), A.cols(), b.rows(), b.cols());

    //p,q must be permutation vectors
    if (p.length() != A.rows())
        throw error::invalid_permvec_length(p.length(), A.rows());

    if (q.length() != A.cols())
        throw error::invalid_permvec_length(q.length(), A.cols());

    Integer N = A.rows();

    value_code vc_A = A.get_value_code();
    value_code vc_B = b.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_A, vc_B);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
    value_code vcr  = matrix_traits::real_value_type(vc);

    if (vc == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_rev2");

    if (b.rows() == 0 || N == 0)
    {
        ret = details::convert_value(Matrix(b), vc);
        return;
    };

    if (b.structural_nnz() == 0)
    {
        ret = details::convert_value(Matrix(b), vc);
        return;
    };

    if (N == 1)
    {
        if (A == 0)
            throw error::error_singular();

        Matrix one  = ones(1,1, vcr);

        if (trans == trans_type::conj_trans)
        {
            ret = b * div(one, matcl::conj(A));
            return;
        }
        else
        {
            ret = b * div(one, A);
            return;
        };
    }

    if (A.get_struct().is_id())
        return linsolve_id_rev(ret, p,q,b, trans, opts, vc);

    Integer ldiags      = matcl::get_ld(A, 0);
    Integer udiags      = matcl::get_ud(A, 0);

    using pivot_type    = opt::linsolve::pivot_type;
    Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
    pivot_type piv      = (pivot_type)piv_int;
    bool isp            = (piv == pivot_type::partial);

    if (ldiags == 0 && udiags == 0)
        return linsolve_diag_rev(ret, A,p,q, b, trans, opts);

    if (isp && ldiags == 0)
        return linsolve_triang_rev(ret, A,p,q, b,false, trans, opts);

    if (isp && udiags == 0)
        return linsolve_triang_rev(ret, A,p,q, b,true, trans, opts);

    if (is_unitary(A.get_struct()))
        return linsolve_unitary_rev(ret, A,p,q, b, trans, opts);

    //quasi-triangular matrices case is not implemented
    return linsolve_general_rev(ret, A,p,q, b, trans, opts);
};

static linsolve_obj linsolve_balanced_impl(const Matrix& A, const Matrix& Dl, const linsolve_obj& B, const Matrix& Dr,
                                           bool has_A, bool sym)
{
    if (Dl.is_vector() == false)
        throw error::vector_required(Dl.rows(), Dl.cols());
    if (sym == false && Dr.is_vector() == false)
        throw error::vector_required(Dr.rows(), Dr.cols());

    if (Dl.length() != B.rows())
        throw error::invalid_mul(Dl.length(), Dl.length(), B.rows(), B.cols(), trans_type::no_trans,
                                 trans_type::no_trans);
    if (sym == false && Dr.length() != B.cols())
        throw error::invalid_mul(B.rows(), B.cols(), Dr.length(), Dr.length(), trans_type::no_trans,
                                 trans_type::no_trans);
    if (has_A == true && (A.rows() != B.rows() || A.cols() != B.cols()) )
        throw error::invalid_size2(B.rows(), B.cols(), A.rows(), A.cols());

    if (Dl.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_balanced");
    if (sym == false && Dr.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_balanced");
    if (has_A == true && A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_balanced");

    bool isv    = Dl.all_finite() && B.all_finite() && (sym == true || Dr.all_finite())
                && (has_A == false || A.all_finite());

    using data_ptr = linsolve_obj::linsolve_data_ptr;

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(Dl.get_value_code(), Dr.get_value_code());
        vc              = matrix_traits::unify_value_types(vc, B.get_value_code());
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

        if (has_A == true)
            vc          = matrix_traits::unify_value_types(vc, A.get_value_code());

        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(B.rows(), vc, B.get_type() )));
    };

    using vec_tag = details::linsolve_obj_diag::from_vec_inv;

    if (sym == false)
    {
        linsolve_obj lo_Dl(data_ptr(new details::linsolve_obj_diag(Dl,vec_tag())));
        linsolve_obj lo_Dr(data_ptr(new details::linsolve_obj_diag(Dr,vec_tag())));

        if (has_A == false)
            return linsolve_seq(lo_Dl, B, lo_Dr);
        else
            return linsolve_seq(A, lo_Dl, B, lo_Dr);
    }
    else
    {
        linsolve_obj lo_Dl(data_ptr(new details::linsolve_obj_diag(Dl,vec_tag())));

        if (has_A == false)
            return linsolve_seq(lo_Dl, B, lo_Dl);
        else
            return linsolve_seq(A, lo_Dl, B, lo_Dl);
    };
};

};};

namespace matcl
{

linsolve_obj matcl::linsolve_lu(const Matrix& A0, const matcl::options& opts)
{
    Matrix A(A0);

    linsolve_obj ret;
    details::linsolve_lu_impl(ret, A, opts);
    return ret;
};

linsolve_obj matcl::linsolve_lu(Matrix&& A0, const matcl::options& opts)
{
    Matrix A(std::move(A0));

    linsolve_obj ret;
    details::linsolve_lu_impl(ret, A, opts);
    return ret;
};

Matrix matcl::linsolve(const matcl::Matrix& A0, const matcl::Matrix& b0, trans_type trans,
                       const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(b0);

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve(const matcl::Matrix& A0, matcl::Matrix&& b0, trans_type trans,
                       const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(std::move(b0));

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve(matcl::Matrix&& A0, const matcl::Matrix& b0, trans_type trans,
                       const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(b0);

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve(matcl::Matrix&& A0, matcl::Matrix&& b0, trans_type trans,
                       const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(std::move(b0));

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve(const matcl::Matrix& A0, const matcl::permvec& p, 
                const matcl::permvec& q, const matcl::Matrix& b0, trans_type trans,
                const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(b0);

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve(const matcl::Matrix& A0, const matcl::permvec& p, 
                const matcl::permvec& q, matcl::Matrix&& b0, trans_type trans,
                const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(std::move(b0));

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve(matcl::Matrix&& A0, const matcl::permvec& p, 
                const matcl::permvec& q, const matcl::Matrix& b0, trans_type trans,
                const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(b0);

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve(matcl::Matrix&& A0, const matcl::permvec& p, 
                const matcl::permvec& q, matcl::Matrix&& b0, trans_type trans,
                const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(std::move(b0));

    Matrix ret;
    details::linsolve_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev(const matcl::Matrix& A0,const matcl::Matrix& b0, const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(b0);

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev(const matcl::Matrix& A0, matcl::Matrix&& b0, const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(std::move(b0));

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev(matcl::Matrix&& A0,const matcl::Matrix& b0, const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(b0);

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev(matcl::Matrix&& A0, matcl::Matrix&& b0, const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(std::move(b0));

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev(const matcl::Matrix& A0, const matcl::permvec& p, 
                           const matcl::permvec& q, const matcl::Matrix& b0,
                           const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(b0);
    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev(const matcl::Matrix& A0, const matcl::permvec& p, 
                           const matcl::permvec& q, matcl::Matrix&& b0,
                           const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(std::move(b0));
    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev(matcl::Matrix&& A0, const matcl::permvec& p, 
                           const matcl::permvec& q, const matcl::Matrix& b0,
                           const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(b0);
    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev(matcl::Matrix&& A0, const matcl::permvec& p, 
                           const matcl::permvec& q, matcl::Matrix&& b0,
                           const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(std::move(b0));
    Matrix ret;
    details::linsolve_rev_impl(ret, A,p,q,b, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(const matcl::Matrix& A0,const matcl::Matrix& b0, trans_type trans,
                            const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(b0);

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(const matcl::Matrix& A0, matcl::Matrix&& b0, trans_type trans,
                            const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(std::move(b0));

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(matcl::Matrix&& A0,const matcl::Matrix& b0, trans_type trans,
                            const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(b0);

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(matcl::Matrix&& A0, matcl::Matrix&& b0, trans_type trans,
                            const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(std::move(b0));

    permvec p = permvec::identity(A.rows());
    permvec q = permvec::identity(A.cols());

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(const matcl::Matrix& A0, const matcl::permvec& p, 
                            const matcl::permvec& q,
                            const matcl::Matrix& b0, trans_type trans,
                            const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(b0);

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(const matcl::Matrix& A0, const matcl::permvec& p, 
                            const matcl::permvec& q,
                            matcl::Matrix&& b0, trans_type trans, const matcl::options& opts)
{
    Matrix A(A0);
    Matrix b(std::move(b0));

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(matcl::Matrix&& A0, const matcl::permvec& p, 
                            const matcl::permvec& q,
                            const matcl::Matrix& b0, trans_type trans,
                            const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(b0);

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::linsolve_rev2(matcl::Matrix&& A0, const matcl::permvec& p, 
                            const matcl::permvec& q,
                            matcl::Matrix&& b0, trans_type trans, const matcl::options& opts)
{
    Matrix A(std::move(A0));
    Matrix b(std::move(b0));

    Matrix ret;
    details::linsolve_rev2_impl(ret, A,p,q,b,trans, opts);
    return ret;
};

Matrix matcl::inv(const matcl::Matrix& A0, const matcl::options& opts)
{
    return make_linsolve_obj(A0, opts).inv();
};

Matrix matcl::inv(matcl::Matrix&& A0, const matcl::options& opts)
{
    return make_linsolve_obj(std::move(A0), opts).inv();
};

linsolve_obj matcl::linsolve_nan(Integer N, value_code vc)
{
    ti::ti_object tiv;
    switch (vc)
    {
        case value_code::v_integer:
            tiv = ti::ti_object_type<Integer>();
            break;
        case value_code::v_real:
            tiv = ti::ti_object_type<Real>();
            break;
        case value_code::v_float:
            tiv = ti::ti_object_type<Float>();
            break;
        case value_code::v_complex:
            tiv = ti::ti_object_type<Complex>();
            break;
        case value_code::v_float_complex:
            tiv = ti::ti_object_type<Float_complex>();
            break;
        case value_code::v_object:
            throw error::object_value_type_not_allowed("linsolve_object");
        default:
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
    }

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, tiv )));
};

linsolve_obj matcl::linsolve_unitary(const Matrix& U)
{
    if (U.rows() != U.cols())
        throw error::square_matrix_required(U.rows(), U.cols());

    if (U.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_unitary");

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    
    if (U.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(U.rows(), U.get_value_code(), U.get_type() )));
    };

    if (U.get_struct().is_id())
        return linsolve_obj(data_ptr(new details::linsolve_obj_id(U.rows(), U.get_value_code(), U.get_type())));

    return linsolve_obj(data_ptr(new details::linsolve_obj_unitary(U)));
};

linsolve_obj matcl::linsolve_unitary(const unitary_matrix& U)
{
    if (U.get_impl()->is_matrix())
        return linsolve_unitary(U.to_matrix());

    if (U.rows() != U.cols())
        throw error::square_matrix_required(U.rows(), U.cols());

    if (U.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_unitary");

    if (U.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(U.rows(), U.get_value_code(), U.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_unitary_mat(U)));
};

linsolve_obj matcl::linsolve_diag(const Matrix& D, const options& opts)
{
    if (D.rows() != D.cols())
        throw error::square_matrix_required(D.rows(), D.cols());

    if (matcl::is_diag(D) == false)
        throw error::diagonal_matrix_required(matcl::get_ld(D,-1), matcl::get_ud(D,-1));

    if (D.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_diag");

    using data_ptr = linsolve_obj::linsolve_data_ptr;

    if (D.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(D.rows(), D.get_value_code(), D.get_type() )));
    };

    if (D.get_struct().is_id())
        return linsolve_obj(data_ptr(new details::linsolve_obj_id(D.rows(), D.get_value_code(), D.get_type())));

    return linsolve_obj(data_ptr(new details::linsolve_obj_diag(D, opts)));
};

linsolve_obj matcl::linsolve_balanced(const Matrix& Dl, const linsolve_obj& B, const Matrix& Dr)
{
    return linsolve_balanced_impl(0, Dl, B, Dr, false, false);
};

linsolve_obj matcl::linsolve_balanced(const Matrix& A, const Matrix& Dl, const linsolve_obj& B, const Matrix& Dr)
{
    return linsolve_balanced_impl(A, Dl, B, Dr, true, false);
};

linsolve_obj matcl::linsolve_balanced_sym(const Matrix& Dlr, const linsolve_obj& B)
{
    return linsolve_balanced_impl(0, Dlr, B, Dlr, false, true);
};

linsolve_obj matcl::linsolve_balanced_sym(const Matrix& A, const Matrix& Dlr, const linsolve_obj& B)
{
    return linsolve_balanced_impl(A, Dlr, B, Dlr, true, true);
};

linsolve_obj matcl::linsolve_perm(const Matrix& A, const linsolve_obj& Apq, const permvec& p,
                                      const permvec& q)
{
    if (A.rows() != Apq.rows() || A.cols() != Apq.cols())
        throw error::invalid_size2(A.rows(), A.cols(), Apq.rows(), Apq.cols());

    if (p.length() != Apq.rows())
        throw error::invalid_length_of_permutation(p.length(), Apq.rows());
    if (q.length() != Apq.cols())
        throw error::invalid_length_of_permutation(q.length(), Apq.cols());

    if (Apq.all_finite() == false)
        return Apq;

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_perm(A, Apq, p, q, false)));
};

linsolve_obj matcl::linsolve_perm(const linsolve_obj& Apq, const permvec& p, const permvec& q)
{
    if (p.length() != Apq.rows())
        throw error::invalid_length_of_permutation(p.length(), Apq.rows());
    if (q.length() != Apq.cols())
        throw error::invalid_length_of_permutation(q.length(), Apq.cols());

    if (Apq.all_finite() == false)
        return Apq;

    if (p.is_id() == true && q.is_id() == true)
        return Apq;

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_perm(Apq, p, q, false)));
};

linsolve_obj matcl::linsolve_symperm(const Matrix& A, const linsolve_obj& Apq, const permvec& p)
{
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    if (A.rows() != Apq.rows() || A.cols() != Apq.cols())
        throw error::invalid_size2(A.rows(), A.cols(), Apq.rows(), Apq.cols());

    if (p.length() != Apq.rows())
        throw error::invalid_length_of_permutation(p.length(), Apq.rows());

    if (Apq.all_finite() == false)
        return Apq;

    if (p.is_id() == true)
        return Apq;

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_perm(A, Apq, p, p, true)));
};

linsolve_obj matcl::linsolve_symperm(const linsolve_obj& Apq, const permvec& p)
{
    if (Apq.rows() != Apq.cols())
        throw error::square_matrix_required(Apq.rows(), Apq.cols());

    if (p.length() != Apq.rows())
        throw error::invalid_length_of_permutation(p.length(), Apq.rows());

    if (Apq.all_finite() == false)
        return Apq;

    if (p.is_id() == true)
        return Apq;

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_perm(Apq, p, p, true)));
};

linsolve_obj matcl::linsolve_extended_sol(const linsolve_obj& lo, const options& opts)
{
    if (lo.all_finite() == false)
        return lo;

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_extended_sol(lo, opts)));
}

linsolve_obj matcl::linsolve_diag_22(const Matrix& D, const options& opts)
{
    Integer ld  = matcl::get_ld(D,1);
    Integer ud  = matcl::get_ud(D,1);

    if (ld == 0 && ud == 0)
        return linsolve_diag(D, opts);

    if (D.rows() != D.cols())
        throw error::square_matrix_required(D.rows(), D.cols());

    if (ld > 1 || ud > 1)
        throw error::band_matrix_required(matcl::get_ld(D,-1), matcl::get_ud(D,-1), 1, 1);

    if (D.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_diag_22");

    if (D.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(D.rows(), D.get_value_code(), D.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_diag_22(D, opts)));
};

linsolve_obj matcl::linsolve_triang(const Matrix& T, const options& opts)
{
    if (T.rows() != T.cols())
        throw error::square_matrix_required(T.rows(), T.cols());

    Integer ld  = matcl::get_ld(T,0);
    Integer ud  = matcl::get_ud(T,0);

    if (ld > 0 && ud > 0)
        throw error::triangular_matrix_required(matcl::get_ld(T,-1), matcl::get_ud(T,-1));

    if (T.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_triang");

    using data_ptr = linsolve_obj::linsolve_data_ptr;

    if (ld == 0 && ud == 0)
        return linsolve_diag(T, opts);

    permvec id = permvec::identity(T.rows());

    if (T.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(T.rows(), T.get_value_code(), T.get_type() )));
    };

    return linsolve_obj(data_ptr(new details::linsolve_obj_triang(T,id,id, opts)));
}

linsolve_obj matcl::linsolve_triang(const Matrix& T, const permvec& p, const permvec& q, const options& opts)
{
    if (T.rows() != T.cols())
        throw error::square_matrix_required(T.rows(), T.cols());

    if (matcl::get_ld(T,0) > 0 && matcl::get_ud(T,0) > 0)
        throw error::triangular_matrix_required(matcl::get_ld(T,-1), matcl::get_ud(T,-1));

    if (p.length() != T.rows())
        throw error::invalid_length_of_permutation(p.length(), T.rows());
    if (p.length() != T.rows() || q.length() != T.cols())
        throw error::invalid_length_of_permutation(q.length(), T.cols());

    if (T.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_triang");

    if (p.is_id() && q.is_id() && is_diag(T))
        return linsolve_diag(T, opts);

    if (T.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(T.rows(), T.get_value_code(), T.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_triang(T,p,q, opts)));
};

linsolve_obj matcl::linsolve_uhess(const Matrix& T, const options& opts)
{
    Integer ld = matcl::get_ld(T,1);

    if (ld == 0)
        return linsolve_triang(T, opts);

    if (T.rows() != T.cols())
        throw error::square_matrix_required(T.rows(), T.cols());

    if (ld > 1)
        throw error::band_matrix_required(matcl::get_ld(T,-1), matcl::get_ud(T,-1), 1, -1);

    if (T.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_uhess");

    using data_ptr = linsolve_obj::linsolve_data_ptr;

    if (T.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(T.rows(), T.get_value_code(), T.get_type() )));
    };
    
    return linsolve_obj(data_ptr(new details::linsolve_obj_uhess(T, opts)));
};

linsolve_obj matcl::linsolve_lhess(const Matrix& T, const options& opts)
{
    Integer ud = matcl::get_ud(T,1);

    if (ud == 0)
        return linsolve_triang(T, opts);

    if (T.rows() != T.cols())
        throw error::square_matrix_required(T.rows(), T.cols());

    if (ud > 1)
        throw error::band_matrix_required(matcl::get_ld(T,-1), matcl::get_ud(T,-1), -1, 1);

    if (T.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_lhess");

    using data_ptr = linsolve_obj::linsolve_data_ptr;

    if (T.all_finite() == false)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(T.rows(), T.get_value_code(), T.get_type() )));
    };

    permvec id = permvec::identity(T.rows());
    
    return linsolve_obj(data_ptr(new details::linsolve_obj_lhess_dense(T, opts)));
};

linsolve_obj matcl::linsolve_seq(const linsolve_obj& A1, const linsolve_obj& A2)
{
    if (A1.cols() != A2.rows())
        throw error::invalid_mul(A1.rows(), A1.cols(), A2.rows(), A2.cols(), trans_type::no_trans,
                                 trans_type::no_trans);

    if (A1.rows() != A1.cols())
        throw error::square_matrix_required(A1.rows(), A1.cols());
    if (A2.rows() != A2.cols())
        throw error::square_matrix_required(A2.rows(), A2.cols());

    if (A1.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");
    if (A2.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");

    if (A1.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A1.get_type())));
    };

    bool isv            = A1.all_finite() && A2.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(A1.rows(), vc, A1.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_seq_2(A1, A2)));
}

linsolve_obj matcl::linsolve_seq(const Matrix& A, const linsolve_obj& A1, const linsolve_obj& A2)
{
    if (A1.cols() != A2.rows())
        throw error::invalid_mul(A1.rows(), A1.cols(), A2.rows(), A2.cols(), trans_type::no_trans,
                                 trans_type::no_trans);

    if (A.rows() != A1.rows())
        throw error::invalid_size2(A1.rows(), A1.cols(), A.rows(), A1.cols());
    if (A.cols() != A2.cols())
        throw error::invalid_size2(A2.rows(), A2.cols(), A2.rows(), A.cols());

    if (A.rows() != A.cols())
        throw error::square_matrix_required(A1.rows(), A1.cols());
    if (A1.rows() != A1.cols())
        throw error::square_matrix_required(A1.rows(), A1.cols());
    if (A2.rows() != A2.cols())
        throw error::square_matrix_required(A2.rows(), A2.cols());

    if (A1.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");
    if (A2.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");
    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");

    if (A1.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A1.get_type())));
    };

    bool isv            = A1.all_finite() && A2.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(A1.rows(), vc, A1.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_seq_2(A, A1, A2)));
}

linsolve_obj matcl::linsolve_seq(const linsolve_obj& A1, const linsolve_obj& A2, 
                                                 const linsolve_obj& A3)
{
    if (A1.cols() != A2.rows())
        throw error::invalid_mul(A1.rows(), A1.cols(), A2.rows(), A2.cols(), trans_type::no_trans,
                                 trans_type::no_trans);
    if (A2.cols() != A3.rows())
        throw error::invalid_mul(A2.rows(), A2.cols(), A3.rows(), A3.cols(), trans_type::no_trans,
                                 trans_type::no_trans);

    if (A1.rows() != A1.cols())
        throw error::square_matrix_required(A1.rows(), A1.cols());
    if (A2.rows() != A2.cols())
        throw error::square_matrix_required(A2.rows(), A2.cols());
    if (A3.rows() != A3.cols())
        throw error::square_matrix_required(A3.rows(), A3.cols());

    if (A1.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");
    if (A2.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");
    if (A3.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");

    if (A1.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A3.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A1.get_type())));
    };

    bool isv            = A1.all_finite() && A2.all_finite() && A3.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A3.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(A1.rows(), vc, A1.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_seq_3(A1, A2, A3)));
}

linsolve_obj matcl::linsolve_seq(const Matrix& A, const linsolve_obj& A1, const linsolve_obj& A2, 
                                                 const linsolve_obj& A3)
{
    if (A1.cols() != A2.rows())
        throw error::invalid_mul(A1.rows(), A1.cols(), A2.rows(), A2.cols(), trans_type::no_trans,
                                 trans_type::no_trans);
    if (A2.cols() != A3.rows())
        throw error::invalid_mul(A2.rows(), A2.cols(), A3.rows(), A3.cols(), trans_type::no_trans,
                                 trans_type::no_trans);

    if (A.rows() != A1.rows())
        throw error::invalid_size2(A1.rows(), A1.cols(), A.rows(), A1.cols());
    if (A.cols() != A3.cols())
        throw error::invalid_size2(A3.rows(), A3.cols(), A3.rows(), A.cols());

    if (A1.rows() != A1.cols())
        throw error::square_matrix_required(A1.rows(), A1.cols());
    if (A2.rows() != A2.cols())
        throw error::square_matrix_required(A2.rows(), A2.cols());
    if (A3.rows() != A3.cols())
        throw error::square_matrix_required(A3.rows(), A3.cols());

    if (A1.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");
    if (A2.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");
    if (A3.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_seq");

    if (A1.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A3.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A1.get_type())));
    };

    bool isv            = A1.all_finite() && A2.all_finite() && A3.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(A1.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(A2.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A3.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(A1.rows(), vc, A1.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_seq_3(A, A1, A2, A3)));
}

tuple<Matrix,bool,bool> matcl::make_nonsingular(const Matrix& T0, Real tol)
{
    Matrix T(T0);

    if (T.rows() != T.cols())
        throw error::square_matrix_required(T.rows(), T.cols());

    Integer ld  = get_ld(T0, 0);
    Integer ud  = get_ud(T0, 0);

    if (ld > 0 && ud > 0)
        throw error::triangular_matrix_required(get_ld(T0,-1), get_ud(T0,-1));

    tuple<Matrix,bool,bool> ret;
    details::vis_make_nonsingular::make<const Matrix&>(T, ret, tol);
    return ret;
};
tuple<Matrix,bool,bool> matcl::make_nonsingular(Matrix&& T0, Real tol)
{
    Matrix T(std::move(T0));

    if (T.rows() != T.cols())
        throw error::square_matrix_required(T.rows(), T.cols());

    Integer ld  = get_ld(T0, 0);
    Integer ud  = get_ud(T0, 0);

    if (ld > 0 && ud > 0)
        throw error::triangular_matrix_required(get_ld(T0,-1), get_ud(T0,-1));

    tuple<Matrix,bool,bool> ret;
    details::vis_make_nonsingular::make<const Matrix&>(T, ret, tol);
    return ret;
};

};
