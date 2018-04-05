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

#pragma once

#include "matcl-core/matrix/enums.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace raw { namespace details
{

template<class SP1, class SP2>
Integer estim_mult_nnz(const SP1 &A, const SP2 &B, trans_type t_A, trans_type t_B)
{
    Integer M, K, N;

    if (t_B == trans_type::no_trans)
        N   = B.cols();
    else
        N   = B.rows();

    if (t_A == trans_type::no_trans)
    {
        M = A.rows();
        K = A.cols();
    }
    else
    {
        K = A.rows();
        M = A.cols();
    };

    Integer A_nnz = A.nnz();
    Integer B_nnz = B.nnz();    

    double da   = double(A_nnz) / ( double(M) * double(K) );
    double db   = double(B_nnz) / ( double(N) * double(K) );
    double pb0  = mrd::pow_helper<double,double>::eval(1. - db, K);

    if (pb0 == 1.)
        return 0;

    double rb   = double(B_nnz) / double(N);
    rb          = (rb/(1-pb0) < K)? rb/(1-pb0) : K;
    
    double ed   = (1. - matcl::pow(1. - da, rb))*(1-pb0);

    return icast( double(M) * double(N) * ed * 1.05);
};

template<bool Val_real>
static const char* trans_code(trans_type t)
{
    switch(t)
    {
        case trans_type::no_trans:
            return "N";
        case trans_type::trans:
            return "T";
        case trans_type::conj_trans:
            return (Val_real)? "T" : "C";
    };

    return nullptr;
};

template<bool Conj_X>
struct make_conj
{
    template<class T>
    static const T& eval(const T& v)
    {
        return v;
    };
};

template<>
struct make_conj<true>
{
    template<class T>
    static T eval(const T& v)
    {
        return conj(v);
    };
};

// alpha operator to form alpha * op(A) * op(B)
// eval Y = Y + alpha * x
template<class Val_ret>
struct alpha_one
{
    template<class Val>
    Val_ret eval(const Val_ret& Y, const Val& v) const
    {
        return Y + v;
    };

    template<class Val>
    Val_ret eval(const Val& v) const
    {
        return v;
    };

    Val_ret get_alpha() const 
    { 
        return Val_ret(1.0); 
    };
};

template<class Val_ret>
struct alpha_mone
{
    template<class Val>
    Val_ret eval(const Val_ret& Y, const Val& v) const
    {
        return Y - v;
    };

    template<class Val>
    Val_ret eval(const Val& v) const
    {
        return -v;
    };

    Val_ret get_alpha() const 
    { 
        return Val_ret(-1.0); 
    };
};

template<class Val_ret>
struct alpha_val
{
    const Val_ret&  alpha;

    alpha_val(const Val_ret& a) : alpha(a) {};

    template<class Val>
    Val_ret eval(const Val_ret& Y, const Val& v) const
    {
        return Y + alpha*v;
    };

    template<class Val>
    Val_ret eval(const Val& v) const
    {
        return alpha*v;
    };

    Val_ret get_alpha() const 
    { 
        return alpha; 
    };
};

//form C = beta * C
template<class Val>
struct prepare_gemm_C
{
    using Mat   = raw::Matrix<Val,struct_dense>;
    static void eval(const Val& beta, Mat& C, Integer fr, Integer rows, Integer fc, Integer cols)
    {
        Integer ld_C    = C.ld();
        Val* ptr_C      = C.ptr() + fr + fc * ld_C;

        level1::ay_test_mat<true, Val, Val, 0, 0, 0>::eval(ptr_C, ld_C, rows, cols, beta);
    };
};

}}}
