/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matmult/func/raw/raw_shprod.h"
#include "matcl-internals/base/instantiate.h"

#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/lib_functions/func_matrix.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-internals/func/inplace.h"

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/algs/scatter.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-matmult/func/raw/mmul/mmul.h"
#include "matcl-internals/func/raw_manip.h"
#include "matcl-internals/func/raw_func_unary.h"
#include "matcl-internals/func/raw_func_plus.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace raw { namespace details
{

namespace md    = matcl::details;
namespace mrd   = matcl::raw::details;
namespace mdyf  = matcl::dynamic::functions;

static trans_type get_mult_trans_type(bool trans, bool is_herprod, int pos)
{
    if (is_herprod == false)
    {
        if (trans == false)
            return pos == 1 ? trans_type::no_trans : trans_type::trans;
        else
            return pos == 2 ? trans_type::no_trans : trans_type::trans;
    }
    else
    {
        if (trans == false)
            return pos == 1 ? trans_type::no_trans : trans_type::conj_trans;
        else
            return pos == 1 ? trans_type::no_trans : trans_type::conj_trans;
    };
};

//------------------------------------------------------------------------
//                      HELPERS
//------------------------------------------------------------------------
template<class VT, bool is_herprod>
struct aat_helper
{
    static VT eval(const VT& val)
    {
        return val * val;
    };
};

template<class VT>
struct aat_helper<VT,true>
{
    static VT eval(const VT& val)
    {
        return mrd::abs2_helper<VT>::eval(val);
    };
};

template<class VT, bool is_herprod>
struct ctrans_helper
{
    static VT eval(const VT& val)
    {
        return val;
    };
};

template<class VT>
struct ctrans_helper<VT,true>
{
    static VT eval(const VT& val)
    {
        return mrd::conj_helper<VT>::eval(val);
    };
};


template<class VT, bool is_herprod>
struct asat_helper
{
    using VTR   = typename md::real_type<VT>::type;

    static VT eval(const VT& val)
    {
        return val * VTR(2);
    };
};

template<>
struct asat_helper<Object,false>
{
    static Object eval(const Object& val)
    {
        return val + ctrans_helper<Object,false>::eval(val);
    };
};

template<>
struct asat_helper<Integer, true>
{
    static Integer eval(Integer val)
    {
        return val * 2;
    };
};
template<>
struct asat_helper<Real, true>
{
    static Real eval(Real val)
    {
        return val * 2.0;
    };
};
template<>
struct asat_helper<Float, true>
{
    static Float eval(Float val)
    {
        return val * 2.0f;
    };
};
template<>
struct asat_helper<Complex, true>
{
    static Real eval(Complex val)
    {
        return real(val) * 2.0;
    };
};
template<>
struct asat_helper<Float_complex, true>
{
    static Float eval(Float_complex val)
    {
        return real(val) * 2.0f;
    };
};
template<>
struct asat_helper<Object, true>
{
    static Object eval(const Object& val)
    {
        return val + ctrans_helper<Object,true>::eval(val);
    };
};

template<class ret_type, class val_type, class struct_type>
struct symprod_helper2 {};

//------------------------------------------------------------------------
//                      DENSE
//------------------------------------------------------------------------
template<class VT>
struct symprod_dense_lapack
{
    using matrix_type = raw::Matrix<VT, struct_dense>;

    static void eval_symprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {        
        const char* str_trans;
        Integer N, K;

        if (trans == false)
        {
            N           = m.rows();
            K           = m.cols();
            str_trans   = "No";
        }
        else
        {
            N           = m.cols();
            K           = m.rows();
            str_trans   = "Trans";
        };

        matrix_type res(ti::ti_empty(), N,N);

        VT One = 1.;
        VT Zero = 0.;

        lapack::syrk("upper", str_trans, N, K, *md::lap(&One), md::lap(m.ptr()), m.ld(), *md::lap(&Zero),
                    md::lap(res.ptr()), res.ld());

        //lower triangular part
        matcl::raw::inplace::copy_utril_to_ltril<matrix_type>::eval(res, false);

        ret = matcl::Matrix(res,true);
        return;
    };
    static void eval_herprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {        
        using VTR = typename md::real_type<VT>::type;

        const char* str_trans;
        Integer N, K;
        if (trans == false)
        {
            N           = m.rows();
            K           = m.cols();
            str_trans   = "No";
        }
        else
        {
            N           = m.cols();
            K           = m.rows();
            str_trans   = "Trans";
        };

        matrix_type res(ti::ti_empty(), N,N);

        VTR One     = 1.;
        VTR Zero    = 0.;

        lapack::herk("upper", str_trans, N, K, *md::lap(&One), md::lap(m.ptr()), m.ld(), *md::lap(&Zero),
                    md::lap(res.ptr()), res.ld());

        //lower triangular part
        matcl::raw::inplace::copy_utril_to_ltril<matrix_type>::eval(res, true);

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class VT, bool is_herprod>
struct symprod_dense_generic
{
    using matrix_type = raw::Matrix<VT, struct_dense>;

    static void eval_sum(matcl::Matrix& ret, const matrix_type& A)
    {
        using ti_ret_type   = typename ti::get_ti_type<VT>::type;
        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_plus::eval(), 
                                                            ti::get_ti(A), ti::get_ti(A));

        Integer N       = A.rows();

        matrix_type C(ret_ti, N, N);

        Integer A_ld    = A.ld();
        Integer C_ld    = C.ld();

        //block size
        const Integer NB    = optim_params::block_size_trans;

        //blocked version, reduce cache misses in matrix transpose
        //prepare upper triangular part only
        for (Integer ib = 0; ib < N; ib += NB)
        {
            Integer max_jb      = std::min(N, ib + NB);

            for (Integer jb = 0; jb < max_jb; jb += NB)
            {
                Integer ibl     = std::min(N, ib + NB);

                const VT* ptr_A = A.ptr() + ib * A_ld;
                const VT* ptr_At= A.ptr() + ib;
                VT * ptr_C      = C.ptr() + ib * C_ld;

                for (Integer i = ib; i < ibl; ++i)
                {                    
                    for (Integer j = jb; j <= i; ++j)
                    {
                        VT vt       = ctrans_helper<VT,is_herprod>::eval(ptr_At[j*A_ld]);
                        ptr_C[j]    = ptr_A[j] + vt;
                    };

                    ptr_C       += C_ld;
                    ptr_A       += A_ld;
                    ptr_At      += 1;
                };
            };
        };

        //lower triangular part
        matcl::raw::inplace::copy_utril_to_ltril<matrix_type>::eval(C, is_herprod);

        ret = matcl::Matrix(C,true);
        if (is_herprod)
            ret.add_struct(predefined_struct_type::her);
        else
            ret.add_struct(predefined_struct_type::sym);
    };

    static void eval(matcl::Matrix& ret, const matrix_type& A, bool trans)
    {
        using ti_ret_type   = typename ti::get_ti_type<VT>::type;
        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            ti::get_ti(A), ti::get_ti(A));

        static const bool is_compl  = md::is_complex<VT>::value;
        VT Z                        = md::default_value<VT>(ret_ti);

        Integer N, K;
        if (trans == false)
        {
            N = A.rows();
            K = A.cols();
        }
        else
        {
            N = A.cols();
            K = A.rows();
        };

        matrix_type C(ret_ti, N, N);

        Integer A_ld    = A.ld();
        Integer C_ld    = C.ld();

        if (trans == false)
        {            
            VT * ptr_C = C.ptr();

            for (Integer j = 0; j < N; ++j)
            {
                for (Integer i = 0; i <= j; ++i)
                    ptr_C[i] = Z;

                const VT* ptr_A = A.ptr();

                for (Integer l = 0; l < K; ++l)
                {
                    VT tmp = ptr_A[j];

                    if (mrd::is_zero(tmp) == false)
                    {
                        tmp = ctrans_helper<VT,is_herprod>::eval(tmp);

                        for (Integer i = 0; i <= j; ++i)
                            ptr_C[i] = ptr_C[i] + tmp * ptr_A[i];
                    };

                    ptr_A += A_ld;
                };

                ptr_C += C_ld;
            };
        }
        else
        {
            const VT* ptr_Aj    = A.ptr();
            VT * ptr_C          = C.ptr();

            for (Integer j = 0; j < N; ++j)
            {
                const VT* ptr_Ai = A.ptr();

                for (Integer i = 0; i <= j; ++i)
                {
                    VT tmp = Z;
                    for (Integer l = 0; l < K; ++l)
                    {
                        VT tr   = ctrans_helper<VT,is_herprod>::eval(ptr_Ai[l]);
                        tmp     = tmp + tr * ptr_Aj[l];
                    };

                    ptr_C[i]    = tmp;
                    ptr_Ai      += A_ld;
                };

                ptr_Aj  += A_ld;
                ptr_C   += C_ld;
            };
        };

        //clear imaginary part on the main diagonal
        if (is_herprod == true && is_compl == true)
        {
            Integer s           = N;
            VT * ptr_C          = C.ptr();

            for (Integer i = 0; i < s; ++i)
            {
                ptr_C[0]        = real(ptr_C[0]);
                ptr_C           += C_ld + 1;
            };
        };

        //lower triangular part
        matcl::raw::inplace::copy_utril_to_ltril<matrix_type>::eval(C, is_herprod);

        ret = matcl::Matrix(C,true);
        return;
    };
};

template<class ret_type, class Val, bool Is_float>
struct symprod_helper_dense
{
    using matrix_type = raw::Matrix<Val, struct_dense>;

    static void eval_symprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        return symprod_dense_generic<Val, false>::eval(ret,m,trans);
    };
    static void eval_herprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        return symprod_dense_generic<Val, true>::eval(ret,m,trans);
    };
};

template<class ret_type, class Val>
struct symprod_helper_dense<ret_type, Val, true>
{
    using matrix_type = raw::Matrix<Val, struct_dense>;    

    static void eval_symprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        symprod_dense_lapack<Val>::eval_symprod(ret,m,trans);

        trans_type tr_1 = get_mult_trans_type(trans, false, 1);
        trans_type tr_2 = get_mult_trans_type(trans, false, 2);
        bool is_square  = ret.rows() == ret.cols();
        bool is_sq_mat  = m.rows() == m.cols();
        value_code vc   = matrix_traits::value_code<Val>::value;
        struct_flag sf  = predefined_struct_ext::seft_mult_struct(m.get_struct(), tr_1, tr_2,
                                vc, is_sq_mat, is_square);

        ret.add_struct(sf);
        return;
    }
    static void eval_herprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        symprod_dense_lapack<Val>::eval_herprod(ret,m,trans);

        trans_type tr_1 = get_mult_trans_type(trans, true, 1);
        trans_type tr_2 = get_mult_trans_type(trans, true, 2);
        bool is_square  = ret.rows() == ret.cols();
        bool is_sq_mat  = m.rows() == m.cols();
        value_code vc   = matrix_traits::value_code<Val>::value;
        struct_flag sf  = predefined_struct_ext::seft_mult_struct(m.get_struct(), tr_1, tr_2,
                                vc, is_sq_mat, is_square);

        ret.add_struct(sf);
        return;
    };
};

template<class ret_type, class val_type>
struct symprod_helper2<ret_type,val_type,struct_dense>
    :symprod_helper_dense<ret_type,val_type, md::is_float_scalar<val_type>::value>
{
    static void eval_symsum(matcl::Matrix& ret, const matrix_type& m)
    {
        return symprod_dense_generic<val_type, false>::eval_sum(ret,m);
    };
    static void eval_hersum(matcl::Matrix& ret, const matrix_type& m)
    {
        return symprod_dense_generic<val_type, true>::eval_sum(ret,m);
    };
};

//------------------------------------------------------------------------
//                      BAND
//------------------------------------------------------------------------
template<class ret_type, class val_type, bool is_herprod>
struct symprod_helper_band
{
    using matrix_type = raw::Matrix<val_type, struct_banded>;

    static void eval_sum(matcl::Matrix& ret, const matrix_type& A)
    {
        using VT            = val_type;
        using ti_ret_type   = typename ti::get_ti_type<VT>::type;
        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_plus::eval(), 
                                                            ti::get_ti(A), ti::get_ti(A));

        Integer N   = A.rows();

        Integer fd  = A.first_diag();
        Integer ld  = A.last_diag();

        Integer fdr = std::min(fd, -ld);
        Integer ldr = std::max(ld, -fd);

        val_type Z  = md::default_value<val_type>(ret_ti);

        bool has_d0 = A.has_diag(0);

        ret_type C  = has_d0 ? ret_type(ret_ti, N, N, fdr, ldr) 
                             : ret_type(ret_ti, Z, N, N, fdr, ldr);

        Integer A_ld    = A.ld();
        Integer C_ld    = C.ld();

        //upper triangular part
        for (Integer d = 1; d <= ldr; ++d)
        {
            Integer s       = C.diag_length(d);

            if (A.has_diag(d) && A.has_diag(-d))
            {               
                const VT* A_ptr     = A.rep_ptr() + A.first_elem_diag(d);
                const VT* At_ptr    = A.rep_ptr() + A.first_elem_diag(-d);
                VT* C_ptr           = C.rep_ptr() + C.first_elem_diag(d);

                for (Integer i = 0; i < s; ++i)
                {
                    VT at           = ctrans_helper<VT, is_herprod>::eval(At_ptr[0]);
                    C_ptr[0]        = A_ptr[0] + at;
                    C_ptr           += C_ld;
                    A_ptr           += A_ld;
                    At_ptr          += A_ld;
                };
            }
            else if (A.has_diag(d))
            {
                const VT* A_ptr     = A.rep_ptr() + A.first_elem_diag(d);
                VT* C_ptr           = C.rep_ptr() + C.first_elem_diag(d);

                for (Integer i = 0; i < s; ++i)
                {
                    C_ptr[0]        = A_ptr[0];
                    C_ptr           += C_ld;
                    A_ptr           += A_ld;
                };
            }
            else if (A.has_diag(-d))
            {
                const VT* At_ptr    = A.rep_ptr() + A.first_elem_diag(-d);
                VT* C_ptr           = C.rep_ptr() + C.first_elem_diag(d);

                for (Integer i = 0; i < s; ++i)
                {
                    VT at           = ctrans_helper<VT, is_herprod>::eval(At_ptr[0]);
                    C_ptr[0]        = at;
                    C_ptr           += C_ld;
                    At_ptr          += A_ld;
                };
            }
        };

        //diagonal
        if (has_d0)
        {
            Integer s           = C.diag_length(0);

            const VT* A_ptr     = A.rep_ptr() + A.first_elem_diag(0);
            VT* C_ptr           = C.rep_ptr() + C.first_elem_diag(0);

            for (Integer i = 0; i < s; ++i)
            {
                val_type sum    = asat_helper<val_type, is_herprod>::eval(A_ptr[0]);
                C_ptr[0]        = sum;
                C_ptr           += C_ld;
                A_ptr           += A_ld;
            };
        }

        //lower triangular part
        matcl::raw::inplace::copy_utril_to_ltril<matrix_type>::eval(C, is_herprod);

        ret = matcl::Matrix(C,true);
        if (is_herprod)
            ret.add_struct(predefined_struct_type::her);
        else
            ret.add_struct(predefined_struct_type::sym);
    };

    static void eval(matcl::Matrix& ret, const matrix_type& A, bool trans)
    {
        using ti_ret_type  = typename ti::get_ti_type<val_type>::type;
        ti_ret_type ret_ti = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            ti::get_ti(A), ti::get_ti(A));

        Integer r   = A.rows();
        Integer c   = A.cols();

        Integer N;

        static const bool is_compl  = md::is_complex<val_type>::value;
        val_type Z                  = md::default_value<val_type>(ret_ti);

        if (trans == false)
            N   = r;
        else
            N   = c;

        Integer fd  = A.first_diag();
        Integer ld  = A.last_diag();

        Integer l   = std::min(std::max(N-1, 0), ld - fd);
        Integer u   = l;

        ret_type C(ret_ti, N, N, -l, u);

        Integer A_ld    = A.ld();
        Integer C_ld    = C.ld();

        //upper triangular part
        if (trans == false)
        {            
            const val_type* ptr_A = A.rep_ptr();
            val_type * ptr_C      = C.rep_ptr();

            for (Integer j = 0; j < N; ++j)
            {
                Integer ic      = C.first_elem_pos(j);
                Integer A_fc    = A.first_col(j);
                Integer A_lc    = A.last_col(j);
                Integer max_i   = std::min(j, C.last_row(j));

                for (Integer i = C.first_row(j); i <= max_i; ++i, ++ic)
                {				
                    ptr_C[ic]   = Z;

                    Integer ks  = std::max(A.first_col(i), A_fc);
                    Integer ke  = std::min(A.last_col(i), A_lc);

                    Integer ia0 = A.first_elem_pos(ks) - A.first_row(ks) + imult(ks, A_ld);
                    Integer ia  = ia0 + i;
                    Integer ib  = ia0 + j;

                    for (Integer k = ks; k <= ke; ++k, ia += A_ld-1, ib += A_ld-1)
                    {
                        val_type B_tr = ctrans_helper<val_type, is_herprod>::eval(ptr_A[ib]);
                        ptr_C[ic]   = ptr_C[ic] + ptr_A[ia] * B_tr;
                    };
                };

                ptr_C   += C_ld;
            };
        }
        else
        {
            const val_type* ptr_A0  = A.rep_ptr();
            const val_type* ptr_B   = ptr_A0;
            val_type * ptr_C        = C.rep_ptr();

            for (Integer j = 0; j < N; ++j)
            {
                Integer ic      = C.first_elem_pos(j);
                Integer B_fr    = A.first_row(j);
                Integer B_lr    = A.last_row(j);
                Integer ib0     = A.first_elem_pos(j) - A.first_row(j);
                Integer i       = C.first_row(j);
                Integer max_i   = std::min(j, C.last_row(j));

                const val_type* ptr_A = ptr_A0 + i*A_ld;

                for (; i <= max_i; ++i, ++ic)
                {				
                    ptr_C[ic]   = Z;

                    Integer ks  = std::max(A.first_row(i), B_fr);
                    Integer ke  = std::min(A.last_row(i), B_lr);
                    Integer ia  = A.first_elem_pos(i) - A.first_row(i) + ks;
                    Integer ib  = ib0 + ks;

                    for (Integer k = ks; k <= ke; ++k, ++ia, ++ib)
                    {
                        val_type A_tr   = ctrans_helper<val_type, is_herprod>::eval(ptr_A[ia]);
                        ptr_C[ic]       = ptr_C[ic] + A_tr * ptr_B[ib];
                    };

                    ptr_A   += A_ld;
                };

                ptr_B   += A_ld;
                ptr_C   += C_ld;
            };
        };

        //lower triangular part
        matcl::raw::inplace::copy_utril_to_ltril<ret_type>::eval(C, is_herprod);

        //clear imaginary part on the main diagonal
        if (is_herprod == true && is_compl == true && C.has_diag(0) == true)
        {
            Integer s           = C.diag_length(0);
            val_type * ptr_C    = C.rep_ptr() + C.first_elem_diag(0);

            for (Integer i = 0; i < s; ++i)
            {
                ptr_C[0]        = real(ptr_C[0]);
                ptr_C           += C_ld;
            };
        };

        ret = matcl::Matrix(C,true);

        trans_type tr_1 = get_mult_trans_type(trans, is_herprod, 1);
        trans_type tr_2 = get_mult_trans_type(trans, is_herprod, 2);
        bool is_square  = ret.rows() == ret.cols();
        bool is_sq_A    = A.rows() == A.cols();
        value_code vc   = matrix_traits::value_code<val_type>::value;
        struct_flag sf  = predefined_struct_ext::seft_mult_struct(A.get_struct(), tr_1, tr_2,
                                vc, is_sq_A, is_square);

        ret.add_struct(sf);
        return;
    };
};

template<class ret_type, class val_type>
struct symprod_helper2<ret_type, val_type, struct_banded>
{
    using matrix_type = raw::Matrix<val_type, struct_banded>;

    static void eval_symsum(matcl::Matrix& ret, const matrix_type& m)
    {
        return symprod_helper_band<ret_type,val_type, false>::eval_sum(ret,m);
    };
    static void eval_hersum(matcl::Matrix& ret, const matrix_type& m)
    {
        return symprod_helper_band<ret_type,val_type, true>::eval_sum(ret,m);
    };

    static void eval_symprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        return symprod_helper_band<ret_type,val_type, false>::eval(ret,m,trans);
    };
    static void eval_herprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        return symprod_helper_band<ret_type,val_type, true>::eval(ret,m,trans);
    };
};

//------------------------------------------------------------------------
//                      SPARSE
//------------------------------------------------------------------------
template<class ret_type, class val_type, bool is_herprod>
struct symprod_helper_sparse
{
    using matrix_type = raw::Matrix<val_type, struct_sparse>;

    static void eval_sum(matcl::Matrix& ret, const matrix_type& A)
    {
        matrix_type At(A.get_type());

        // hard to find enything better
        if (is_herprod)
            manip_trans_helper<matrix_type>::eval_ctrans(At, A);
        else
            manip_trans_helper<matrix_type>::eval_trans(At, A);

        plus_helper_mat_mat_inpl<matrix_type, matrix_type>::eval(ret, A, At);

        if (is_herprod)
            ret.add_struct(predefined_struct_type::her);
        else
            ret.add_struct(predefined_struct_type::sym);
    };

    static void eval(matcl::Matrix& ret, const matrix_type& A, bool trans)
    {
        if (trans == false)
        {
            if (is_herprod == false)
            {
                matcl::Matrix At    = matcl::trans(matcl::Matrix(A,false));
                matrix_type A2      = At.get_impl<matrix_type>();
                return eval(ret,A2,true);
            }
            else
            {
                matcl::Matrix At    = matcl::ctrans(matcl::Matrix(A,false));
                matrix_type A2      = At.get_impl<matrix_type>();
                return eval(ret, A2,true);
            };
        };
        
        const sparse_ccs<val_type>& Ad = A.rep();
        Integer N = A.cols();
        Integer K = A.rows();

        const Integer * A_c     = Ad.ptr_c();
        const Integer * A_r     = Ad.ptr_r();
        const val_type* A_x     = Ad.ptr_x();

        using ti_ret_type   = typename ti::get_ti_type<val_type>::type;
        using workspace     = md::workspace2<val_type>;

        ti_ret_type ret_ti  = ti::get_return_ti<ti_ret_type>(mdyf::op_mul::eval(), 
                                                            ti::get_ti(A), ti::get_ti(A));

        matcl::algorithm::scatter sc_i      = matcl::algorithm::scatter::get(K,2*N);
        workspace sc_v                      = workspace(ret_ti,K);

        matcl::pod_workspace<Integer> work_ind(N,0);

        //estimate nnz
        Integer nz = 0;
        for (Integer j = 0; j < N; ++j)
        {
            Integer i1  = A_c[j];
            Integer i2  = A_c[j+1];

            if (i2 == i1)
                continue;

            auto mark       = sc_i.next_mark();

            for (Integer i = i1; i < i2; ++i)
            {
                Integer k   = A_r[i];
                sc_i[k]     = mark;
            };

            for (Integer i = 0; i < j; ++i)
            {
                Integer j1  = A_c[i];
                Integer j2  = A_c[i+1];

                if (j2 == j1)
                    continue;

                for (Integer k = j1; k < j2; ++k)
                {
                    Integer l = A_r[k];

                    if (sc_i[l] == mark)
                    {
                        ++work_ind[i];
                        ++work_ind[j];
                        nz += 2;
                        break;
                    };
                };    
            };

            ++nz;
            ++work_ind[j];
        };

        //make product
        matrix_type C(ret_ti, N, N, nz);

        sparse_ccs<val_type>& Cd     = C.rep();

        Integer * C_c       = Cd.ptr_c();
        Integer * C_r       = Cd.ptr_r();
        val_type* C_x       = Cd.ptr_x();
        val_type Z          = md::default_value<val_type>(ret_ti);

        C_c[0] = 0;
        for (Integer j = 1; j <= N; ++j)
        {
            C_c[j]          = C_c[j-1] + work_ind[j-1];
            work_ind[j-1]   = 0;
        };

        for (Integer j = 0; j < N; ++j)
        {
            Integer i1  = A_c[j];
            Integer i2  = A_c[j+1];

            if (i2 == i1)
                continue;

            auto mark       = sc_i.next_mark();

            for (Integer i = i1; i < i2; ++i)
            {
                Integer k   = A_r[i];
                sc_i[k]     = mark;
                sc_v[k]     = ctrans_helper<val_type, is_herprod>::eval(A_x[i]);
            };

            for (Integer i = 0; i < j; ++i)
            {
                Integer j1  = A_c[i];
                Integer j2  = A_c[i+1];

                if (j2 == j1)
                    continue;

                val_type val = Z;
                bool any     = false;

                for (Integer k = j1; k < j2; ++k)
                {
                    Integer l = A_r[k];

                    if (sc_i[l] == mark)
                    {
                        any     = true;
                        val     = val + sc_v[l] * A_x[k];
                    };
                };

                if (any == false)
                    continue;
                
                {
                    Integer& pos        = work_ind[j];
                    C_r[C_c[j] + pos]   = i;
                    C_x[C_c[j] + pos]   = ctrans_helper<val_type, is_herprod>::eval(val);
                    ++pos;
                };
                {
                    Integer& pos        = work_ind[i];
                    C_r[C_c[i] + pos]   = j;
                    C_x[C_c[i] + pos]   = val;
                    ++pos;
                };
            };

            //diagonal element
            val_type val = Z;

            for (Integer k = i1; k < i2; ++k)
                val   = val + aat_helper<val_type, is_herprod>::eval(A_x[k]);
            
            if (is_herprod == false)
            {
                Integer& pos        = work_ind[j];
                C_r[C_c[j] + pos]   = j;
                C_x[C_c[j] + pos]   = val;
                ++pos;
            }
            else
            {
                Integer& pos        = work_ind[j];
                C_r[C_c[j] + pos]   = j;
                C_x[C_c[j] + pos]   = real(val);
                ++pos;
            }
        };

        ret = matcl::Matrix(C,true);

        trans_type tr_1 = get_mult_trans_type(trans, is_herprod, 1);
        trans_type tr_2 = get_mult_trans_type(trans, is_herprod, 2);
        bool is_square  = ret.rows() == ret.cols();
        bool is_sq_mat  = A.rows() == A.cols();
        value_code vc   = matrix_traits::value_code<val_type>::value;
        struct_flag sf  = predefined_struct_ext::seft_mult_struct(A.get_struct(), tr_1, tr_2, 
                                vc, is_sq_mat, is_square);

        ret.add_struct(sf);
        return;
    };
};

template<class ret_type, class val_type>
struct symprod_helper2<ret_type, val_type, struct_sparse>
{
    using matrix_type = raw::Matrix<val_type, struct_sparse>;

    static void eval_symsum(matcl::Matrix& ret, const matrix_type& m)
    {
        return symprod_helper_sparse<ret_type,val_type,false>::eval_sum(ret,m);
    };
    static void eval_hersum(matcl::Matrix& ret, const matrix_type& m)
    {
        return symprod_helper_sparse<ret_type,val_type,true>::eval_sum(ret,m);
    };

    static void eval_symprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        return symprod_helper_sparse<ret_type,val_type, false>::eval(ret,m,trans);
    };
    static void eval_herprod(matcl::Matrix& ret, const matrix_type& m, bool trans)
    {
        return symprod_helper_sparse<ret_type,val_type, true>::eval(ret,m,trans);
    };
};

//------------------------------------------------------------------------
//                      GENERAL
//------------------------------------------------------------------------
template<class ret_type, class val_type, class str_type>
struct symprod_helper_diag
{
    using matrix_type = raw::Matrix<val_type, str_type>;

    static void eval(matcl::Matrix& ret, const matrix_type& mat, bool trans, bool is_herprod)
    {        
        Integer M   = mat.rows();
        Integer N   = mat.cols();
        Integer K   = (trans == false)? M : N;

        matcl::Matrix diag = matcl::Matrix(mat.get_diag(), false);

        if (is_herprod == true)
            diag    = matcl::abs2(std::move(diag));
        else
            diag    = matcl::mul(diag, std::move(diag));

        diag.resize(K, 1);

        ret = matcl::bdiag(diag);
        
        trans_type tr_1 = get_mult_trans_type(trans, is_herprod, 1);
        trans_type tr_2 = get_mult_trans_type(trans, is_herprod, 2);
        bool is_square  = ret.rows() == ret.cols();
        bool is_sq_mat  = mat.rows() == mat.cols();
        value_code vc   = matrix_traits::value_code<val_type>::value;
        struct_flag sf  = predefined_struct_ext::seft_mult_struct(mat.get_struct(), tr_1, tr_2,
                                vc, is_sq_mat, is_square);
        ret.add_struct(sf);
    };

    static void eval_symsum(matcl::Matrix& ret, const matrix_type& mat)
    {
        matcl::Matrix diag  = matcl::Matrix(mat.get_diag(), false);
        diag                = matcl::mul(std::move(diag), val_type(2.0));

        ret = matcl::bdiag(diag);
    }
    static void eval_hersum(matcl::Matrix& ret, const matrix_type& mat)
    {
        matcl::Matrix diag  = matcl::Matrix(mat.get_diag(), false);
        diag                = matcl::mul(matcl::real(std::move(diag)), val_type(2.0));

        ret = matcl::bdiag(diag);
    }
};

template<class MP>
void symprod_helper<MP>::eval_symprod(matcl::Matrix& ret, const MP& m, bool trans)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;
    using ret_type  = ret_type_symprod;

    Integer M       = m.rows();
    Integer N       = m.cols();
    Integer K       = (trans == false)? M : N;

    if (m.get_struct().is_id())
    {
        ret = matcl::Matrix(m,false);
        return;
    }

    if (m.get_struct().is_diag())
        return symprod_helper_diag<ret_type, val_type, str_type>::eval(ret,m,trans,false);

    if (is_unitary(m.get_struct()))
    {
        if(std::is_same<val_type,Real>::value == true
            ||std::is_same<val_type,Float>::value == true)
        {
            using val_type_re   = typename md::real_type<val_type>::type;
            auto res            = md::unit_matrix<val_type_re,struct_sparse>::eval(m.get_type(),K);
            ret                 = matcl::Matrix(res,false);
            return;
        };
    }

    if (K == 0 || m.nnz() == 0)
    {
        using val_type_re   = typename md::real_type<val_type>::type;
        auto res            = md::zero_matrix<val_type_re,struct_sparse>::eval(m.get_type(),K,K);
        ret                 = matcl::Matrix(res,false);
        return;
    };

    if (matcl::is_diag(matcl::Matrix(m,false)))
        return symprod_helper_diag<ret_type, val_type, str_type>::eval(ret,m,trans,false);

    return symprod_helper2<ret_type, val_type, str_type>::eval_symprod(ret,m,trans);
};

template<class MP>
void symprod_helper<MP>::eval_herprod(matcl::Matrix& ret, const MP& m, bool trans)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;
    using ret_type  = ret_type_symprod;

    Integer M       = m.rows();
    Integer N       = m.cols();
    Integer K       = (trans == false)? M : N;

    if (m.get_struct().is_id())      
    {
        ret = matcl::Matrix(m,false);
        return;
    }
    if (m.get_struct().is_diag())
        return symprod_helper_diag<ret_type, val_type, str_type>::eval(ret,m,trans,true);

    if (is_unitary(m.get_struct()))
    {
        if(std::is_same<val_type,Real>::value == true || std::is_same<val_type,Float>::value == true
            || md::is_complex<val_type>::value == true)
        {
            using val_type_re   = typename md::real_type<val_type>::type;
            auto res            = md::unit_matrix<val_type_re,struct_sparse>::eval(m.get_type(),K);
            ret                 = matcl::Matrix(res,false);
            return;
        };
    }

    if (K == 0 || m.nnz() == 0)
    {
        using val_type_re   = typename md::real_type<val_type>::type;
        auto res            = md::zero_matrix<val_type_re,struct_sparse>::eval(m.get_type(),K,K);
        ret                 = matcl::Matrix(res,false);
        return;
    };

    if (matcl::is_diag(matcl::Matrix(m,false)))
        return symprod_helper_diag<ret_type, val_type, str_type>::eval(ret,m,trans,true);

    return symprod_helper2<ret_type, val_type, str_type>::eval_herprod(ret,m,trans);
};

template<class MP>
void symprod_helper<MP>::eval_symsum(matcl::Matrix& ret, const MP& m)
{
    if (m.rows() != m.cols())
        throw error::square_matrix_required(m.rows(), m.cols());

    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;
    using ret_type  = ret_type_symprod;
    using VR        = typename md::real_type<val_type>::type;

    if (m.get_struct().is_diag())
        return symprod_helper_diag<ret_type, val_type, str_type>::eval_symsum(ret,m);

    if (m.nnz() == 0)
    {
        ret                 = matcl::Matrix(m,false);
        return;
    };

    bool is_real    = md::is_complex<val_type>::value == false;
    bool is_sym     = m.get_struct().is_symmetric(true, is_real);
    bool is_her     = m.get_struct().is_hermitian(true, is_real);

    if (is_sym == true)
    {
        mult_helper_mat_scal<MP, val_type>::eval(ret, m, val_type(2.0), trans_type::no_trans, trans_type::no_trans);
        return;
    };
    if (is_her == true)
    {
        scalfunc_real_helper<MP>::eval_real(ret, m);
        ret = mmul(std::move(ret), VR(2.0));
        return;
    };

    if (matcl::is_diag(matcl::Matrix(m,false)))
        return symprod_helper_diag<ret_type, val_type, str_type>::eval_symsum(ret,m);

    return symprod_helper2<ret_type, val_type, str_type>::eval_symsum(ret,m);
};

template<class MP>
void symprod_helper<MP>::eval_hersum(matcl::Matrix& ret, const MP& m)
{
    if (m.rows() != m.cols())
        throw error::square_matrix_required(m.rows(), m.cols());

    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;
    using ret_type  = ret_type_symprod;
    using VR        = typename md::real_type<val_type>::type;

    if (m.get_struct().is_diag())
        return symprod_helper_diag<ret_type, val_type, str_type>::eval_hersum(ret,m);

    if (m.nnz() == 0)
    {
        ret                 = matcl::Matrix(m,false);
        return;
    };

    bool is_real    = md::is_complex<val_type>::value == false;
    bool is_her     = m.get_struct().is_hermitian(true, is_real);
    bool is_sym     = m.get_struct().is_symmetric(true, is_real);

    if (is_her == true)
    {
        mult_helper_mat_scal<MP, val_type>::eval(ret, m, val_type(2.0), trans_type::no_trans, trans_type::no_trans);
        return;
    };
    if (is_sym == true)
    {
        scalfunc_real_helper<MP>::eval_real(ret, m);
        ret = mmul(std::move(ret), VR(2.0));
        return;
    };

    if (matcl::is_diag(matcl::Matrix(m,false)))
        return symprod_helper_diag<ret_type, val_type, str_type>::eval_hersum(ret,m);

    return symprod_helper2<ret_type, val_type, str_type>::eval_hersum(ret,m);
};

};};};

#pragma warning( pop )

MACRO_INSTANTIATE_G_1(matcl::raw::details::symprod_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::symprod_helper)
