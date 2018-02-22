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

#include "matcl-matmult/func/raw/kron.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-matrep/matrix/struct_flag_ext.h"

namespace matcl { namespace raw { namespace details
{

namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;
namespace md    = matcl::details;
namespace mdyf  = matcl::dynamic::functions;

static struct_flag kron_struct(struct_flag A, struct_flag B, bool is_re_1, bool is_re_2, 
                               bool is_square_1, bool is_square_2, bool is_square_ret)
{
    struct_flag sf = predefined_struct_ext::kron_struct(A, B, is_re_1, is_re_2, is_square_1, 
                                                        is_square_2, is_square_ret);
    return sf;
};

template<class V, bool Is_float = md::is_float_scalar<V>::value>
struct kron_row_col
{
    using mat_type = Matrix<V,struct_dense>;

    static void eval(matcl::Matrix& ret, const mat_type& A, const mat_type& B)
    {
        Integer As  = A.size();
        Integer Bs  = B.size();
        Integer lda = (A.rows() == 1)? A.ld() : 1;
        Integer ldb = (B.rows() == 1)? B.ld() : 1;

        using ret_ti_type   = ti::ti_type<V>;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        mat_type res(ret_ti,Bs,As);
        const V* ptr_A  = A.ptr();
        V* ptr_res      = res.ptr();
        Integer res_ld  = res.ld();

        for (Integer i = 0; i < As; ++i, ptr_A += lda)
        {
            const V* ptr_B  = B.ptr();
            const V& aij    = ptr_A[0];     

            if (mrd::is_zero(aij) == false)
            {
                for (Integer Bi = 0; Bi < Bs; ++Bi)
                {
                    ptr_res[Bi] = aij * ptr_B[Bi];
                };

                ptr_B += ldb;
            }
            else
            {
                for (Integer Bi = 0; Bi < Bs; ++Bi)
                {
                    ptr_res[Bi] = aij;
                };
            };

            ptr_res += res_ld;
        };

        ret = matcl::Matrix(res,true);

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        struct_flag sf  = predefined_struct_ext
                            ::mult_struct(A.get_struct(), B.get_struct(), trans_type::no_trans,
                                trans_type::no_trans, is_real_matrix(A), is_real_matrix(B), 
                                is_sq_A, is_sq_B, ret.is_square());

        ret.set_struct(sf);
        return;
    }
};

template<class Val>
struct kron_row_col<Val, true>
{
    using mat_type = Matrix<Val,struct_dense>;

    static void eval(matcl::Matrix& ret, const mat_type& A, const mat_type& B)
    {
        Integer As = A.size();
        Integer Bs = B.size();

        mat_type res(A.get_type(),0., Bs, As);

        Val alph = 1.;

        lapack::geru(Bs, As, *md::lap(&alph), md::lap(B.ptr()), 1, md::lap(A.ptr()), A.ld(), 
                    md::lap(res.ptr()), Bs);

        ret = matcl::Matrix(res,true);

        bool is_sq_A    = A.rows() == A.cols();
        bool is_sq_B    = B.rows() == B.cols();
        struct_flag sf  = predefined_struct_ext
                            ::mult_struct(A.get_struct(), B.get_struct(), trans_type::no_trans,
                                trans_type::no_trans, is_real_matrix(A), is_real_matrix(B), 
                                is_sq_A, is_sq_B, ret.is_square());

        ret.set_struct(sf);

        return;
    }
};

template<class ret_type,class M1,class M2,class str1, class str2>
struct eval_kron
{};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_dense,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using val_1     = typename M1::value_type;
        using val_2     = typename M2::value_type;
        using val_ret   = typename ret_type::value_type;

        const val_1* ptr_A = A.ptr();
        const val_2* ptr_B = B.ptr();

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols();
        Integer Br = B.rows(), Bc = B.cols();

        if (Ar == 1 && Ac == 1&& Br == 1 && Bc == 1)
        {
            val_ret val = ptr_A[0] * ptr_B[0];
            ret         = val;
            return;
        };

        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer resr = imult_c(Ar,Br);
        Integer resc = imult_c(Ac,Bc);

        if (resr == 0 || resc == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,resr, resc),false);
            return;
        };	

        // A is row vector, B is column vector.
        if (Ar == 1 && Bc == 1)
        {
            return kron_row_col<val_ret>::eval(ret, converter<ret_type,M1>::eval(A),
                                                    converter<ret_type,M2>::eval(B));
        }

        // B is row vector, A is column vector.
        if (Br == 1 && Ac == 1)
        {
            return kron_row_col<val_ret>::eval(ret, converter<ret_type,M2>::eval(B),
                                                    converter<ret_type,M1>::eval(A));
        }

        ret_type res(ret_ti,resr, resc);
        val_ret* ptr_res    = res.ptr();
        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer res_ld      = res.ld();

        for (Integer Aj = 0; Aj < Ac; ++Aj)
        {
            ptr_B = B.ptr();
            for (Integer Bj = 0; Bj < Bc; ++Bj)
            {
                for (Integer Ai = 0, pos_ret = 0; Ai < Ar; ++Ai)
                {
                    const val_1& aij = ptr_A[Ai];

                    if (mrd::is_zero(aij) == false)
                    {
                        for (Integer Bi = 0; Bi < Br; ++Bi, ++pos_ret)
                        {
                            ptr_res[pos_ret] = aij * ptr_B[Bi];
                        };
                    }
                    else
                    {
                        for (Integer Bi = 0; Bi < Br; ++Bi, ++pos_ret)
                        {
                            ptr_res[pos_ret] = val_ret(aij);
                        };
                    };
                };
                ptr_B   += B_ld;
                ptr_res += res_ld;
            };

            ptr_A += A_ld;
        };

        ret = matcl::Matrix(res,true);
        
        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_dense,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        using val_ret       = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        const val_type_1* ptr_A = A.ptr();

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cn == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cn),false);
            return;
        };
     
        const sparse_ccs<val_type_2>& Bd	= B.rep();
        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const val_type_2* Bd_x		= Bd.ptr_x();

        ret_type C(ret_ti,Cr, Cc, Cn);
        sparse_ccs<val_ret>& d = C.rep(); 
     
        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        val_ret* d_x    = d.ptr_x();

        Integer A_ld    = A.ld();

        for (Integer aj = 0, cj = 0, k = 0; aj < Ac; ++aj)
        {
            for (Integer bj = 0; bj < Bc; ++bj, ++cj)
            {
                d_c[cj + 1] = d_c[cj];
                for (Integer ka = 0, arBr = 0; ka < Ar; ++ka, arBr += Br)
                {
                    const val_type_1& aijx = ptr_A[ka];

                    if (mrd::is_zero(aijx) == false)
                    {
                        d_c[cj + 1] += Bd_c[bj + 1] - Bd_c[bj];
                        for (Integer kb = Bd_c[bj]; kb < Bd_c[bj + 1]; ++kb, ++k)
                        {
                            d_r[k] = arBr + Bd_r[kb];
                            d_x[k] = aijx * Bd_x[kb];
                        };
                    };
                };
            };
            ptr_A += A_ld;
        };

        d.add_memory(-1);
        ret = matcl::Matrix(C,true);

        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_dense,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        using VTR   = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cn == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cn),false);
            return;
        };
     
        ret_type C(ret_ti,Cr, Cc, Cn);
        sparse_ccs<VTR>& d  = C.rep(); 
        const VT1* ptr_A    = A.ptr();

        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        VTR* d_x        = d.ptr_x();
        
        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();
        Integer fdb         = B.first_diag();
        Integer ldb         = B.last_diag();
        Integer k           = 0;

        if (fdb == ldb)
        {
            const VT2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fdb);
            Integer sb          = B.diag_length(fdb);
            Integer c0b         = B.first_col_on_diag(fdb);
            Integer r0b         = B.first_row_on_diag(fdb);

            for (Integer aj = 0, cj = 0; aj < Ac; ++aj)
            {
                const VT2* ptr_Bl   = ptr_B;
                Integer rbj         = r0b;

                for (Integer bj = 0; bj < c0b; ++bj, ++cj)
                    d_c[cj]     = k;

                for (Integer bj = c0b; bj < c0b + sb; ++bj, ++cj, ptr_Bl += B_ld, ++rbj)
                {
                    d_c[cj]     = k;

                    const VT2& tmp = ptr_Bl[0];

                    if (mrd::is_zero(tmp) == true)
                        continue;

                    for (Integer ka = 0, arBr = 0; ka < Ar; ++ka, arBr += Br)
                    {
                        const VT1& aijx = ptr_A[ka];

                        if (mrd::is_zero(aijx) == false)
                        {
                            d_r[k] = arBr + rbj;
                            d_x[k] = aijx * tmp;
                            ++k;
                        };
                    };
                }

                for (Integer bj = c0b + sb; bj < Bc; ++bj, ++cj)
                    d_c[cj]     = k;

                ptr_A += A_ld;
            };
        }
        else
        {
            for (Integer aj = 0, cj = 0; aj < Ac; ++aj)
            {
                const VT2* ptr_B = B.rep_ptr();

                for (Integer bj = 0; bj < Bc; ++bj, ++cj, ptr_B += B_ld)
                {
                    d_c[cj]     = k;

                    for (Integer ka = 0, arBr = 0; ka < Ar; ++ka, arBr += Br)
                    {
                        const VT1& aijx = ptr_A[ka];

                        if (mrd::is_zero(aijx) == false)
                        {
                            Integer first_row   = B.first_row(bj);
                            Integer last_row    = B.last_row(bj);
                            Integer pos_B       = B.first_elem_pos(bj);

                            for (Integer kb = first_row; kb <= last_row; ++kb, ++pos_B)
                            {
                                const VT2& tmp = ptr_B[pos_B];

                                if (mrd::is_zero(tmp) == true)
                                    continue;

                                d_r[k] = arBr + kb;
                                d_x[k] = aijx * tmp;
                                ++k;
                            };
                        };
                    };
                };

                ptr_A += A_ld;
            };            
        };

        d_c[Cc] = k;

        d.add_memory(-1);
        ret = matcl::Matrix(C,true);

        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_sparse,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        using val_ret       = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cn == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cn),false);
            return;
        };
     
        const sparse_ccs<val_type_1>& Ad	= A.rep();
        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const val_type_1* Ad_x		= Ad.ptr_x();

        ret_type C(ret_ti,Cr, Cc, Cn); 

        sparse_ccs<val_ret>& d           = C.rep();
        const val_type_2* ptr_B     = B.ptr();
        Integer B_ld                = B.ld();

        Integer* d_c = d.ptr_c();
        Integer* d_r = d.ptr_r();
        val_ret* d_x = d.ptr_x();
     
        for (Integer aj = 0, cj = 0, k = 0; aj < Ac; ++aj)
        {
            ptr_B = B.ptr();
            for (Integer bj = 0; bj < Bc; ++bj, ++cj)
            {
                d_c[cj + 1] = d_c[cj];
                for (Integer ka = Ad_c[aj]; ka < Ad_c[aj + 1]; ++ka)
                {					
                    const val_type_1& aijx = Ad_x[ka];			
                    
                    if (mrd::is_zero(aijx) == true)
                    {
                        continue;
                    };

                    Integer arBr = imult(Ad_r[ka], Br);
                    Integer nz = 0;
                    for (Integer kb = 0; kb < Br; ++kb)
                    {
                        const val_type_2& tmp = ptr_B[kb]; 
                        if (mrd::is_zero(tmp) == false)
                        {
                            d_r[k] = arBr + kb;
                            d_x[k] = aijx * tmp;
                            ++nz;
                            ++k;
                        };
                    };
                    d_c[cj + 1] += nz;
                };

                ptr_B += B_ld;
            };
        };

        d.add_memory(-1);
        ret = matcl::Matrix(C,true);

        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_sparse,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using val_type_1    = typename M1::value_type;
        using val_type_2    = typename M2::value_type;
        using val_ret       = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cn == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cn),false);
            return;
        };
     
        const sparse_ccs<val_type_1>& Ad	= A.rep();
        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const val_type_1* Ad_x		= Ad.ptr_x();

        const sparse_ccs<val_type_2>& Bd = B.rep();
        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const val_type_2* Bd_x		= Bd.ptr_x();

        ret_type C(ret_ti,Cr, Cc, Cn); 
        sparse_ccs<val_ret>& d = C.rep();
     
        Integer* d_c = d.ptr_c();
        Integer* d_r = d.ptr_r();
        val_ret* d_x = d.ptr_x();

        for (Integer aj = 0, cj = 0, k = 0; aj < Ac; ++aj)
        {
            for (Integer bj = 0; bj < Bc; ++bj, ++cj)
            {
                d_c[cj + 1] = d_c[cj];
                for (Integer ka = Ad_c[aj]; ka < Ad_c[aj + 1]; ++ka)
                {					
                    const val_type_1& aijx = Ad_x[ka];
                    if (mrd::is_zero(aijx) == true)
                    {
                        continue;
                    };

                    Integer arBr = imult(Ad_r[ka],Br);
                    d_c[cj + 1] += Bd_c[bj + 1] - Bd_c[bj];
                    for (Integer kb = Bd_c[bj]; kb < Bd_c[bj + 1]; ++kb, ++k)
                    {
                        d_r[k] = arBr + Bd_r[kb];
                        d_x[k] = aijx * Bd_x[kb];
                    };
                };
            };
        };

        d.add_memory(-1);
        ret = matcl::Matrix(C,true);

        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_sparse,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        using VTR   = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cn == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cn),false);
            return;
        };
     
        const sparse_ccs<VT1>& Ad   = A.rep();
        const Integer* Ad_c         = Ad.ptr_c();
        const Integer* Ad_r         = Ad.ptr_r();
        const VT1* Ad_x		        = Ad.ptr_x();

        ret_type C(ret_ti,Cr, Cc, Cn);

        sparse_ccs<VTR>& d  = C.rep();     
        Integer* d_c        = d.ptr_c();
        Integer* d_r        = d.ptr_r();
        VTR* d_x            = d.ptr_x();
        
        Integer B_ld        = B.ld();
        Integer fdb         = B.first_diag();
        Integer ldb         = B.last_diag();

        Integer k           = 0;

        if (fdb == ldb)
        {
            const VT2* ptr_B = B.rep_ptr() + B.first_elem_diag(fdb);
            Integer sb      = B.diag_length(fdb);
            Integer r0b     = B.first_row_on_diag(fdb);
            Integer c0b     = B.first_col_on_diag(fdb);

            for (Integer aj = 0, cj = 0; aj < Ac; ++aj)
            {
                const VT2* ptr_Bl   = ptr_B;
                Integer rbj         = r0b;

                for (Integer bj = 0; bj < c0b; ++bj, ++cj)
                    d_c[cj] = k;

                for (Integer bj = c0b; bj < c0b + sb; ++bj, ++cj, ++rbj, ptr_Bl += B_ld)
                {
                    d_c[cj] = k;

                    const VT2& tmp = ptr_Bl[0];

                    if (mrd::is_zero(tmp) == true)
                        continue;

                    for (Integer ka = Ad_c[aj]; ka < Ad_c[aj + 1]; ++ka)
                    {					
                        const VT1& aijx = Ad_x[ka];

                        if (mrd::is_zero(aijx) == true)
                            continue;

                        Integer arBr = imult(Ad_r[ka],Br);

                        d_r[k]  = arBr + rbj;
                        d_x[k]  = aijx * tmp;
                        ++k;
                    };
                }

                for (Integer bj = c0b + sb; bj < Bc; ++bj, ++cj)
                    d_c[cj]     = k;
            }
        }
        else
        {
            for (Integer aj = 0, cj = 0; aj < Ac; ++aj)
            {
                const VT2* ptr_B    = B.rep_ptr();

                for (Integer bj = 0; bj < Bc; ++bj, ++cj, ptr_B += B_ld)
                {
                    d_c[cj]     = k;

                    for (Integer ka = Ad_c[aj]; ka < Ad_c[aj + 1]; ++ka)
                    {					
                        const VT1& aijx = Ad_x[ka];
                    
                        if (mrd::is_zero(aijx) == true)
                            continue;

                        Integer arBr = imult(Ad_r[ka],Br);
                        
                        Integer first_row   = B.first_row(bj);
                        Integer last_row    = B.last_row(bj);
                        Integer pos_B       = B.first_elem_pos(bj);

                        for (Integer kb = first_row; kb <= last_row; ++kb, ++pos_B)
                        {
                            const VT2& tmp = ptr_B[pos_B];

                            if (mrd::is_zero(tmp) == true)
                                continue;

                            d_r[k] = arBr + kb;
                            d_x[k] = aijx * tmp;
                            ++k;
                        };
                    };
                };
            };
        };

        d_c[Cc] = k;
        d.add_memory(-1);
        ret = matcl::Matrix(C,true);
        
        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_banded,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        using VTR   = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cr == 0 || Cc == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cc),false);
            return;
        };
     
        ret_type C(ret_ti,Cr, Cc, Cn); 

        sparse_ccs<VTR>& d  = C.rep();
        Integer nz          = 0;

        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        VTR* d_x        = d.ptr_x();
        
        const VT2* ptr_B = B.ptr();

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer fda     = A.first_diag();
        Integer lda     = A.last_diag();

        if (fda == lda)
        {	 
            const VT1* ptr_A = A.rep_ptr() + A.first_elem_diag(fda);
            Integer sa  = A.diag_length(fda);
            Integer c0a = A.first_col_on_diag(fda);

            Integer cj  = 0;

            for (Integer aj = 0; aj <c0a; ++aj)
            {
                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                    d_c[cj] = nz;
            };

            for (Integer aj = c0a, arBr = 0; aj < c0a + sa; ++aj, ptr_A += A_ld, arBr += Br)
            {
                const VT1& aijx = ptr_A[0];
      
                if (mrd::is_zero(aijx) == true)
                {
                    for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                        d_c[cj] = nz;

                    continue;
                };

                ptr_B   = B.ptr();

                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                {
                    d_c[cj] = nz;

                    for (Integer kb = 0; kb < Br; ++kb)
                    {
                        const VT2& tmp = ptr_B[kb]; 
                
                        if (!mrd::is_zero(tmp))
                        {
                            d_r[nz] = arBr + kb;
                            d_x[nz] = aijx * tmp;
                            ++nz;
                        };
                    };

                    ptr_B   += B_ld;
                }
            };

            for (Integer aj = c0a + sa; aj < Ac; ++aj)
            {
                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                    d_c[cj] = nz;
            };

            d_c[Cc] = nz;
        }
        else
        {
            const VT1* ptr_A = A.rep_ptr();

            for (Integer aj = 0, cj = 0; aj < Ac; ++aj, ptr_A += A_ld)
            {
                ptr_B = B.ptr();

                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                {
                    d_c[cj] = nz;

                    Integer fr      = A.first_row(aj);
                    Integer lr      = A.last_row(aj);
                    Integer pos_A   = A.first_elem_pos(aj);
                    Integer arBr    = imult(fr,Br);

                    for (Integer ka = fr; ka <= lr; ++ka, ++pos_A, arBr += Br)
                    {					
                        const VT1& aijx = ptr_A[pos_A];
                     
                        if (mrd::is_zero(aijx) == true)
                            continue;

                        for (Integer kb = 0; kb < Br; ++kb)
                        {
                            const VT2& tmp = ptr_B[kb]; 

                            if (!mrd::is_zero(tmp))
                            {
                                d_r[nz] = arBr + kb;
                                d_x[nz] = aijx * tmp;
                                ++nz;
                            };
                        };
                    };

                    ptr_B += B_ld;
                };
            };
        };

        d_c[Cc] = nz;

        d.add_memory(-1);
        ret = matcl::Matrix(C,true);

        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_banded,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        using VTR   = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cn == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cn),false);
            return;
        };

        const sparse_ccs<VT2>& Bd   = B.rep();
        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const VT2* Bd_x		        = Bd.ptr_x();

        ret_type C(ret_ti,Cr, Cc, Cn); 

        sparse_ccs<VTR>& d = C.rep();

        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        VTR* d_x        = d.ptr_x();

        Integer A_ld    = A.ld();
        Integer fda     = A.first_diag();
        Integer lda     = A.last_diag();
        Integer k       = 0;

        if (fda == lda)
        {
            const VT1* ptr_A    = A.rep_ptr() + A.first_elem_diag(fda); 
            Integer sa          = A.diag_length(fda);
            Integer c0a         = A.first_col_on_diag(fda);

            Integer cj  = 0;

            for (Integer aj = 0; aj < c0a; ++ aj)
            {
                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                    d_c[cj]     = k;
            };

            for (Integer aj = c0a, arBr = 0; aj < c0a + sa; ++aj, ptr_A += A_ld, arBr += Br)
            {
                const VT1& aijx = ptr_A[0];

                if (mrd::is_zero(aijx) == true)
                {
                    for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                        d_c[cj] = k;

                    continue;
                };

                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                {
                    d_c[cj]     = k;

                    for (Integer kb = Bd_c[bj]; kb < Bd_c[bj + 1]; ++kb)
                    {
                        d_r[k] = arBr + Bd_r[kb];
                        d_x[k] = aijx * Bd_x[kb];
                        ++k;
                    };
                }
            };
            for (Integer aj = c0a + sa; aj < Ac; ++ aj)
            {
                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                    d_c[cj]     = k;
            };
        }
        else
        {
            const VT1* ptr_A    = A.rep_ptr(); 

            for (Integer aj = 0, cj = 0; aj < Ac; ++aj, ptr_A += A_ld)
            {
                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                {
                    d_c[cj]         = k;

                    Integer fr      = A.first_row(aj);
                    Integer lr      = A.last_row(aj);
                    Integer pos_A   = A.first_elem_pos(aj);
                    Integer arBr    = imult(fr,Br);

                    for (Integer ka = fr; ka <= lr; ++ka, ++pos_A, arBr += Br)
                    {					
                        const VT1& aijx = ptr_A[pos_A];

                        if (mrd::is_zero(aijx) == true)
                            continue;

                        for (Integer kb = Bd_c[bj]; kb < Bd_c[bj + 1]; ++kb)
                        {
                            d_r[k] = arBr + Bd_r[kb];
                            d_x[k] = aijx * Bd_x[kb];
                            ++k;
                        };
                    };
                };
            };
        };

        d_c[Cc] = k;
        d.add_memory(-1);
        ret = matcl::Matrix(C,true);

        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class ret_type,class M1,class M2>
struct eval_kron<ret_type,M1,M2,struct_banded,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M1& A, const M2& B)
    {
        using VT1   = typename M1::value_type;
        using VT2   = typename M2::value_type;
        using VTR   = typename ret_type::value_type;

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::op_mul::eval(),
                                                                 ti::get_ti(A),ti::get_ti(B));

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();
        Integer Br = B.rows(), Bc = B.cols(), Bn = B.nnz();
     
        error::check_size(Ar,Br);
        error::check_size(Ac,Bc);

        Integer Cr = imult_c(Ar,Br);
        Integer Cc = imult_c(Ac,Bc);
        Integer Cn = imult_c(An,Bn);

        if (Cn == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,Cr,Cn),false);
            return;
        };
     
        ret_type C(ret_ti,Cr, Cc, Cn); 
        sparse_ccs<VTR>& d = C.rep();

        Integer* d_c    = d.ptr_c();
        Integer* d_r    = d.ptr_r();
        VTR* d_x        = d.ptr_x();
     
        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();
        Integer fda     = A.first_diag();
        Integer lda     = A.last_diag();
        Integer fdb     = B.first_diag();
        Integer ldb     = B.last_diag();

        Integer k       = 0;

        if (fda == lda && fdb == ldb)
        {
            const VT1* ptr_A = A.rep_ptr()  + A.first_elem_diag(fda);
            const VT2* ptr_B = B.rep_ptr()  + B.first_elem_diag(fdb);

            Integer c0a     = A.first_col_on_diag(fda);
            Integer sa      = A.diag_length(fda);

            Integer r0b     = B.first_row_on_diag(fdb);
            Integer c0b     = B.first_col_on_diag(fdb);
            Integer sb      = B.diag_length(fdb);

            Integer cj      = 0;

            for (Integer aj = 0; aj < c0a; ++ aj)
            {
                for (Integer j = 0; j < Bc; ++j, ++cj)
                    d_c[cj] = k;
            };
                        
            Integer arBr    = 0;

            for (Integer aj = c0a; aj < c0a + sa; ++aj, ptr_A += A_ld, arBr += Br)
            {
                const VT1& aijx = ptr_A[0];
            
                if (mrd::is_zero(aijx) == true)
                {
                    for (Integer j = 0; j < Bc; ++j, ++cj)
                        d_c[cj] = k;

                    continue;
                };

                const VT2* ptr_Bl   = ptr_B;
                Integer rbj         = r0b;

                for (Integer bj = 0; bj < c0b; ++bj, ++cj)
                    d_c[cj] = k;

                for (Integer bj = c0b; bj < c0b + sb; ++bj)
                {
                    d_c[cj] = k;

                    const VT2& tmp = ptr_Bl[0];

                    if (mrd::is_zero(tmp) == false)
                    {
                        d_r[k] = arBr + rbj;
                        d_x[k] = aijx * tmp;
                    
                        ++k;
                    };

                    cj      += 1;
                    ptr_Bl  += B_ld;
                    rbj     += 1;
                };

                for (Integer bj = c0b + sb; bj < Bc; ++bj, ++cj)
                    d_c[cj] = k;
            };
            
            for (Integer aj = c0a + sa; aj < Ac; ++ aj)
            {
                for (Integer j = 0; j < Bc; ++j, ++cj)
                    d_c[cj] = k;
            };            
        }
        else if (fdb == ldb)
        {
            const VT1* ptr_A = A.rep_ptr();
            const VT2* ptr_B = B.rep_ptr() + B.first_elem_diag(fdb);

            Integer r0b     = B.first_row_on_diag(fdb);
            Integer c0b     = B.first_col_on_diag(fdb);
            Integer sb      = B.diag_length(fdb);

            Integer cj      = 0;

            for (Integer aj = 0; aj < Ac; ++aj, ptr_A += A_ld)
            {
                const VT2* ptr_Bl   = ptr_B;

                for (Integer bj = 0; bj < c0b; ++bj, ++cj)
                    d_c[cj] = k;

                Integer rbj         = r0b;

                for (Integer bj = c0b; bj < c0b + sb; ++bj, ++cj, ptr_Bl += B_ld, ++rbj)
                {
                    d_c[cj] = k;

                    const VT2& tmp = ptr_Bl[0];

                    if (mrd::is_zero(tmp) == true)
                        continue;

                    Integer fr_A    = A.first_row(aj);
                    Integer lr_A    = A.last_row(aj);
                    Integer pos_A   = A.first_elem_pos(aj);
                    Integer arBr    = imult(fr_A,Br) + rbj;

                    for (Integer ka = fr_A; ka <= lr_A; ++ka, ++pos_A, arBr += Br)
                    {					
                        const VT1& aijx = ptr_A[pos_A];

                        if (mrd::is_zero(aijx) == true)
                            continue;

                        d_r[k]  = arBr;
                        d_x[k]  = aijx * tmp;
                        ++k;
                    };
                };
                
                for (Integer bj = c0b + sb; bj < Bc; ++bj, ++cj)
                    d_c[cj] = k;
            };
        }
        else if (fda == lda)
        {
            const VT1* ptr_A = A.rep_ptr()  + A.first_elem_diag(fda);
            const VT2* ptr_B = B.rep_ptr();

            Integer c0a     = A.first_col_on_diag(fda);
            Integer sa      = A.diag_length(fda);

            Integer arBr    = 0;
            Integer cj      = 0;

            for (Integer aj = 0; aj < c0a; ++ aj)
            {
                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                    d_c[cj] = k;
            };

            for (Integer aj = c0a; aj < c0a + sa; ++aj, ptr_A += A_ld, arBr += Br)
            {
                const VT1& aijx = ptr_A[0];
            
                if (mrd::is_zero(aijx) == true)
                {
                    for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                        d_c[cj] = k;

                    continue;
                };

                const VT2* ptr_Bl       = ptr_B;

                for (Integer bj = 0; bj < Bc; ++bj, ++cj, ptr_Bl += B_ld)
                {
                    d_c[cj]             = k;
                
                    Integer first_row   = B.first_row(bj);
                    Integer last_row    = B.last_row(bj);
                    Integer pos_B       = B.first_elem_pos(bj);

                    for (Integer kb = first_row; kb <= last_row; ++kb, ++pos_B)
                    {
                        const VT2& tmp  = ptr_Bl[pos_B];

                        if (mrd::is_zero(tmp) == true)
                            continue;

                        d_r[k]  = arBr + kb;
                        d_x[k]  = aijx * tmp;
                        ++k;
                    };
                };
            };

            for (Integer aj = c0a + sa; aj < Ac; ++ aj)
            {
                for (Integer bj = 0; bj < Bc; ++bj, ++cj)
                    d_c[cj]     = k;
            };
        }
        else
        {
            const VT1* ptr_A = A.rep_ptr();
            const VT2* ptr_B = B.rep_ptr();

            for (Integer aj = 0, cj = 0; aj < Ac; ++aj, ptr_A += A_ld)
            {
                ptr_B       = B.rep_ptr();
                
                for (Integer bj = 0; bj < Bc; ++bj, ++cj, ptr_B += B_ld)
                {
                    d_c[cj] = k;

                    Integer fr_A    = A.first_row(aj);
                    Integer lr_A    = A.last_row(aj);
                    Integer pos_A   = A.first_elem_pos(aj);
                    Integer arBr    = imult(fr_A,Br);

                    for (Integer ka = fr_A; ka <= lr_A; ++ka, ++pos_A, arBr += Br)
                    {					
                        const VT1& aijx = ptr_A[pos_A];
              
                        if (mrd::is_zero(aijx) == true)
                            continue;
                    
                        Integer first_row   = B.first_row(bj);
                        Integer last_row    = B.last_row(bj);
                        Integer pos_B       = B.first_elem_pos(bj);

                        for (Integer kb = first_row; kb <= last_row; ++kb, ++pos_B)
                        {
                            const VT2& tmp = ptr_B[pos_B];

                            if (mrd::is_zero(tmp) == true)
                                continue;

                            d_r[k] = arBr + kb;
                            d_x[k] = aijx * tmp;
                            ++k;
                        };
                    };
                };
            };
        };

        d_c[Cc] = k;
        d.add_memory(-1);

        ret = matcl::Matrix(C,true);

        ret.set_struct(kron_struct(A.get_struct(), B.get_struct(), is_real_matrix(A), 
                                   is_real_matrix(B), Ar == Ac, Br == Bc, ret.is_square()));
        return;
    };
};

template<class M1, class M2>
void kron_helper<M1,M2>::eval(matcl::Matrix& ret, const M1& A, const M2& B)
{
    using str_type_1 = typename M1::struct_type;
    using str_type_2 = typename M2::struct_type;

    return eval_kron<ret_type,M1,M2,str_type_1,str_type_2>::eval(ret,A,B);
}

}}}

MACRO_INSTANTIATE_BIN_ALL(matcl::raw::details::kron_helper)
