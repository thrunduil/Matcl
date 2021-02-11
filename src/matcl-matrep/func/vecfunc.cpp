/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/lib_functions/vecfunc.h"
#include "matcl-matrep/details/extract_type_switch.h"
#include "matcl-matrep/func/raw/raw_vecfunc.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/lib_functions/eval_functors.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-internals/base/utils.h"

namespace matcl { namespace details
{

namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template<class T>
struct sum_helper
{
    static T eval(const T& A, int dim);
};

template<class T>
struct cumsum_helper
{
    static T eval(const T& A, int dim);
};

template<class T>
struct prod_helper
{
    static T eval(const T& A, int dim);
};

template<class T>
struct cumprod_helper
{
    static T eval(const T& A, int dim);
};

template<class T>
struct sumsq_helper
{
    static T eval(const T& A, int dim);
};

template<class T>
struct min_helper
{
    using TR = typename md::real_type<T>::type;

    static T eval(const T& A, int dim);
    static TR eval_abs(const T& A, int dim);
};

template<class T>
struct max_helper
{
    using TR = typename md::real_type<T>::type;

    static T eval(const T& A, int dim);
    static TR eval_abs(const T& A, int dim);
};

template<class T>
struct mean_helper
{
    using R = typename md::unify_types<T,Float>::type;

    static R eval(const T& A, int dim);
};	

template<class T>
struct std_helper
{
    using R = typename md::unify_types<T,Float>::type;
    static R eval(const T& A, int dim);
};	

template<class T>
struct all_helper
{
    static bool eval(const T& A, int dim);
    static bool eval(const T& A,int dim, const test_function& t);
};	

template<class T>
struct any_helper
{
    static bool eval(const T& A, int dim);
    static bool eval(const T& A,int dim, const test_function& t);
};

template<class T>
T sum_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return A;
};

template<class T>
T cumsum_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return A;
};

template<class T>
T prod_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return A;
};

template<class T>
T cumprod_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return A;
};

template<class T>
T sumsq_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return mrd::mul_helper<T,T>::eval(A,A);
};

template<class T>
T min_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return A;
};

template<class T>
typename min_helper<T>::TR
min_helper<T>::eval_abs(const T& A, int dim)
{
    error::check_dim(dim);
    return mrd::abs_helper<T>::eval(A);
};

template<class T>
T max_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return A;
};

template<class T>
typename max_helper<T>::TR
max_helper<T>::eval_abs(const T& A, int dim)
{
    error::check_dim(dim);
    return mrd::abs_helper<T>::eval(A);
};

template<class T>
typename mean_helper<T>::R mean_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return A;
};

template<class T>
typename std_helper<T>::R std_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return md::default_value<T>(ti::get_ti(A));
};

template<class T>
bool all_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return mrd::cast_bool_helper<T>::eval(A);
};

template<class T>
bool all_helper<T>::eval(const T& A,int dim, const test_function& t)
{
    error::check_dim(dim);
    return t.eval(A);
};

template<class T>
bool any_helper<T>::eval(const T& A, int dim)
{
    error::check_dim(dim);
    return mrd::cast_bool_helper<T>::eval(A);
};

template<class T>
bool any_helper<T>::eval(const T& A,int dim, const test_function& t)
{
    error::check_dim(dim);
    return t.eval(A);
};

struct eval_nnz : extract_type_switch<void, eval_nnz, true>
{
    template<class T>
    static void eval(const Matrix&, const T& A, Matrix& ret, Integer dim)
    {
        return mrd::vec_manip_helper<T>::eval_nnz(ret, A, dim);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& A, Matrix& ret, Integer dim)
    {
        error::check_dim(dim);
        ret = mrd::is_zero(A)? 0 : 1;
    };
};

struct eval_nnz_vec : extract_type_switch<void, eval_nnz_vec, true>
{
    template<class T>
    static void eval(const Matrix&, const T& A, Integer& ret)
    {
        return mrd::vec_manip_helper<T>::eval_nnz_vec(ret, A);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& A, Integer& ret)
    {
        ret = mrd::is_zero(A)? 0 : 1;
    };
};

struct eval_all : public extract_type_switch<void,eval_all,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        mrd::vec_manip_helper<T>::eval_all(ret,mat,arg1);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = Matrix(all_helper<T>::eval(mat,arg1));
    };
};

struct eval_all_vec : public extract_type_switch<void,eval_all_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, bool& ret)
    {
        mrd::vec_manip_helper<T>::eval_all_vec(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, bool& ret)
    {
        ret = all_helper<T>::eval(mat,1);
    };
};

struct eval_all_t : public extract_type_switch<void,eval_all_t,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret, const test_function& t,Integer dim)
    {
        mrd::vec_manip_helper<T>::eval_all(ret, mat,t, dim);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret, const test_function& t,Integer dim)
    {
        ret = Matrix(all_helper<T>::eval(mat,dim,t));
    };
};

struct eval_all_t_vec : public extract_type_switch<void,eval_all_t_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, bool& ret, const test_function& t)
    {
        mrd::vec_manip_helper<T>::eval_all_vec(ret, mat,t);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, bool& ret, const test_function& t)
    {
        ret = all_helper<T>::eval(mat,1,t);
    };
};

struct eval_sum : public extract_type_switch<void,eval_sum,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_sum(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& dim)
    {
        ret = sum_helper<T>::eval(mat,dim);
    };
};

struct eval_sum_vec : public extract_type_switch<void,eval_sum_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_sum_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = sum_helper<T>::eval(mat,1);
    };
};

struct eval_prod : public extract_type_switch<void,eval_prod,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_prod(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = prod_helper<T>::eval(mat,arg1);
    };
};

struct eval_prod_vec : public extract_type_switch<void,eval_prod_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_prod_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = prod_helper<T>::eval(mat,1);
    };
};

struct eval_cumsum : public extract_type_switch<void,eval_cumsum,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_cumsum(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = cumsum_helper<T>::eval(mat,arg1);
    };
};

struct eval_cumprod : public extract_type_switch<void,eval_cumprod,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_cumprod(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = cumprod_helper<T>::eval(mat,arg1);
    };
};

struct eval_sumsq : public extract_type_switch<void,eval_sumsq,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_sumsq(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = sumsq_helper<T>::eval(mat,arg1);
    };
};

struct eval_sumsq_vec : public extract_type_switch<void,eval_sumsq_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_sumsq_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = sumsq_helper<T>::eval(mat,1);
    };
};

struct eval_min : public extract_type_switch<void,eval_min,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_min(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = min_helper<T>::eval(mat,arg1);
    };
};

struct eval_min_vec : public extract_type_switch<void,eval_min_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_min_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = min_helper<T>::eval(mat,1);
    };
};

struct eval_max : public extract_type_switch<void,eval_max,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_max(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = max_helper<T>::eval(mat,arg1);
    };
};

struct eval_max_vec : public extract_type_switch<void,eval_max_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_max_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = max_helper<T>::eval(mat,1);
    };
};

struct eval_min_abs : public extract_type_switch<void,eval_min_abs,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_min_abs(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = min_helper<T>::eval_abs(mat,arg1);
    };
};

struct eval_min_abs_vec : public extract_type_switch<void,eval_min_abs_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_min_abs_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = min_helper<T>::eval_abs(mat,1);
    };
};

struct eval_max_abs : public extract_type_switch<void,eval_max_abs,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_max_abs(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = max_helper<T>::eval_abs(mat,arg1);
    };
};

struct eval_max_abs_vec : public extract_type_switch<void,eval_max_abs_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_max_abs_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = max_helper<T>::eval_abs(mat,1);
    };
};

struct eval_mean : public extract_type_switch<void,eval_mean,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& mat, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_mean(ret,mat,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = mean_helper<T>::eval(mat,arg1);
    };
};

struct eval_mean_vec : public extract_type_switch<void,eval_mean_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret)
    {
        return mrd::vec_manip_helper<T>::eval_mean_vec(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret)
    {
        ret = mean_helper<T>::eval(mat,1);
    };
};

struct eval_std : public extract_type_switch<void,eval_std,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& mat, Matrix& ret, const S& dim, bool unbiased)
    {
        return mrd::vec_manip_helper<T>::eval_std(ret, mat, dim,unbiased);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1, bool)
    {
        ret = std_helper<T>::eval(mat,arg1);
    };
};

struct eval_std_vec : public extract_type_switch<void,eval_std_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, Matrix& ret, bool unbiased)
    {
        return mrd::vec_manip_helper<T>::eval_std_vec(ret, mat, unbiased);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, bool)
    {
        ret = std_helper<T>::eval(mat,1);
    };
};

struct eval_any : public extract_type_switch<void,eval_any,true> 
{
    template<class T, class S>
    static void eval(const Matrix& , const T& A, Matrix& ret, const S& dim)
    {
        return mrd::vec_manip_helper<T>::eval_any(ret,A,dim);
    };

    template<class T, class S>
    static void eval_scalar(const Matrix& , const T& mat, Matrix& ret, const S& arg1)
    {
        ret = Matrix(any_helper<T>::eval(mat,arg1));
    };
};

struct eval_any_vec : public extract_type_switch<void,eval_any_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, bool& ret)
    {
        return mrd::vec_manip_helper<T>::eval_any_vec(ret,A);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, bool& ret)
    {
        ret = any_helper<T>::eval(mat,1);
    };
};

struct eval_min2 : public extract_type_switch<void,eval_min2,true> 
{
    template<class T, class S1>
    static void eval(const Matrix& , const T& A, matcl::mat_tup_2& ret, const S1& dim)
    {
        Matrix x, i;
        mrd::vec_manip_helper<T>::eval_min2(x,i,A,dim);
        ret = mat_tup_2(x,i);
    };

    template<class T, class S1>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_tup_2& ret, const S1& d)
    {
        error::check_dim(d);
        ret = mat_tup_2(mat,1);
    };
};

struct eval_min2_vec : public extract_type_switch<void,eval_min2_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, matcl::mat_int2& ret)
    {
        Matrix x;
        Integer i, j;
        mrd::vec_manip_helper<T>::eval_min2_vec(x,i,j,A);
        ret = mat_int2(x,i,j);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_int2& ret)
    {
        ret = mat_int2(mat,1,1);
    };
};

struct eval_max2 : public extract_type_switch<void,eval_max2,true> 
{
    template<class T, class S1>
    static void eval(const Matrix& , const T& mat, matcl::mat_tup_2& ret, const S1& dim)
    {
        Matrix x, i;
        mrd::vec_manip_helper<T>::eval_max2(x,i,mat,dim);
        ret = mat_tup_2(x,i);
    };

    template<class T, class S1>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_tup_2& ret, const S1& d)
    {
        error::check_dim(d);
        ret = mat_tup_2(mat,1);
    };
};

struct eval_max2_vec : public extract_type_switch<void,eval_max2_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, matcl::mat_int2& ret)
    {
        Matrix x;
        Integer i, j;
        mrd::vec_manip_helper<T>::eval_max2_vec(x,i,j,A);
        ret = mat_int2(x,i,j);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_int2& ret)
    {
        ret = mat_int2(mat,1,1);
    };
};

struct eval_min_abs2 : public extract_type_switch<void,eval_min_abs2,true> 
{
    template<class T, class S1>
    static void eval(const Matrix& , const T& mat, matcl::mat_tup_2& ret, const S1& dim)
    {
        Matrix x, i;
        mrd::vec_manip_helper<T>::eval_min_abs2(x,i,mat,dim);
        ret = mat_tup_2(x,i);
    };

    template<class T, class S1>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_tup_2& ret, const S1& d)
    {
        error::check_dim(d);
        ret = mat_tup_2(mrd::abs_helper<T>::eval(mat), 1);
    };
};

struct eval_min_abs2_vec : public extract_type_switch<void,eval_min_abs2_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, matcl::mat_int2& ret)
    {
        Matrix x;
        Integer i, j;
        mrd::vec_manip_helper<T>::eval_min_abs2_vec(x,i,j,A);
        ret = mat_int2(x,i,j);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_int2& ret)
    {
        ret = mat_int2(mrd::abs_helper<T>::eval(mat),1,1);
    };
};

struct eval_max_abs2 : public extract_type_switch<void,eval_max_abs2,true> 
{
    template<class T, class S1>
    static void eval(const Matrix& , const T& mat, matcl::mat_tup_2& ret, const S1& dim)
    {
        Matrix x, i;
        mrd::vec_manip_helper<T>::eval_max_abs2(x,i,mat,dim);
        ret = mat_tup_2(x,i);
    };

    template<class T, class S1>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_tup_2& ret, const S1& d)
    {
        error::check_dim(d);
        ret = mat_tup_2(mrd::abs_helper<T>::eval(mat), 1);
    };
};

struct eval_max_abs2_vec : public extract_type_switch<void,eval_max_abs2_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& A, matcl::mat_int2& ret)
    {
        Matrix x;
        Integer i, j;
        mrd::vec_manip_helper<T>::eval_max_abs2_vec(x,i,j,A);
        ret = mat_int2(x,i,j);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::mat_int2& ret)
    {
        ret = mat_int2(mrd::abs_helper<T>::eval(mat),1,1);
    };
};

struct eval_any_t : public extract_type_switch<void,eval_any_t,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, matcl::Matrix& ret, const test_function& t,Integer dim)
    {
        return mrd::vec_manip_helper<T>::eval_any(ret,mat,t,dim);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, matcl::Matrix& ret, const test_function& t,Integer arg2)
    {
        ret = Matrix(any_helper<T>::eval(mat,arg2,t));
    };
};

struct eval_any_t_vec : public extract_type_switch<void,eval_any_t_vec,true> 
{
    template<class T>
    static void eval(const Matrix& , const T& mat, bool& ret, const test_function& t)
    {
        return mrd::vec_manip_helper<T>::eval_any_vec(ret,mat,t);
    };

    template<class T>
    static void eval_scalar(const Matrix& , const T& mat, bool& ret, const test_function& t)
    {
        ret = any_helper<T>::eval(mat,1,t);
    };
};

}}

namespace matcl
{

Matrix matcl::nnz(const Matrix& A, Integer dim)
{
    matcl::Matrix ret;
    details::eval_nnz::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::all(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_all::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::sum(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_sum::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::cumsum(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_cumsum::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::prod(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_prod::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::cumprod(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_cumprod::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::min_d(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_min::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::max_d(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_max::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::min_abs_d(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_min_abs::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::max_abs_d(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_max_abs::make<const Matrix&>(A,ret,dim);
    return ret;
};

mat_tup_2 matcl::min2(const Matrix& A,int dim)
{
    mat_tup_2 ret;
    details::eval_min2::make<const Matrix&>(A,ret,dim);
    return ret;
};

mat_tup_2 matcl::max2(const Matrix& A,int dim)
{
    mat_tup_2 ret;
    details::eval_max2::make<const Matrix&>(A,ret,dim);
    return ret;
};

mat_tup_2 matcl::min_abs2(const Matrix& A,int dim)
{
    mat_tup_2 ret;
    details::eval_min_abs2::make<const Matrix&>(A,ret,dim);
    return ret;
};

mat_tup_2 matcl::max_abs2(const Matrix& A,int dim)
{
    mat_tup_2 ret;
    details::eval_max_abs2::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::mean(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_mean::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::std(const Matrix& A,int dim,bool unbiased)
{
    Matrix ret;
    details::eval_std::make<const Matrix&>(A,ret,dim,unbiased);
    return ret;
};

Matrix matcl::any(const Matrix& A,int dim)
{
    Matrix ret;
    details::eval_any::make<const Matrix&>(A,ret,dim);
    return ret;
};

Matrix matcl::any(const Matrix& A,const test_function& t,int dim)
{
    Matrix ret;
    details::eval_any_t::make<const Matrix&>(A,ret,t,dim);
    return ret;
};

Matrix matcl::all(const Matrix& A,const test_function& t,int dim)
{
    Matrix ret;
    details::eval_all_t::make<const Matrix&>(A,ret,t,dim);
    return ret;
};

Integer matcl::nnz_vec(const Matrix& A)
{
    Integer ret;
    details::eval_nnz_vec::make<const Matrix&>(A,ret);
    return ret;
};

bool matcl::all_vec(const Matrix& v)
{
    bool ret;
    details::eval_all_vec::make<const Matrix&>(v,ret);
    return ret;
};

bool matcl::all_vec(const Matrix& v,const test_function& t)
{
    bool ret;
    details::eval_all_t_vec::make<const Matrix&>(v,ret,t);
    return ret;
};

bool matcl::any_vec(const Matrix& v)
{
    bool ret;
    details::eval_any_vec::make<const Matrix&>(v,ret);
    return ret;
};

bool matcl::any_vec(const Matrix& v,const test_function& t)
{
    bool ret;
    details::eval_any_t_vec::make<const Matrix&>(v,ret,t);
    return ret;
};

Matrix matcl::sum_vec(const Matrix& v)
{
    Matrix ret;
    details::eval_sum_vec::make<const Matrix&>(v,ret);
    return ret;
};

Matrix matcl::prod_vec(const Matrix& v)
{
    Matrix ret;
    details::eval_prod_vec::make<const Matrix&>(v,ret);
    return ret;
};

Matrix matcl::mean_vec(const Matrix& v)
{
    Matrix ret;
    details::eval_mean_vec::make<const Matrix&>(v,ret);
    return ret;
};

Matrix matcl::std_vec(const Matrix& v, bool unbiased)
{
    Matrix ret;
    details::eval_std_vec::make<const Matrix&>(v,ret,unbiased);
    return ret;
};

Matrix matcl::min_vec(const Matrix& v)
{
    Matrix ret;
    details::eval_min_vec::make<const Matrix&>(v,ret);
    return ret;
};

Matrix matcl::max_vec(const Matrix& v)
{
    Matrix ret;
    details::eval_max_vec::make<const Matrix&>(v,ret);
    return ret;
};

Matrix matcl::min_abs_vec(const Matrix& v)
{
    Matrix ret;
    details::eval_min_abs_vec::make<const Matrix&>(v,ret);
    return ret;
};

Matrix matcl::max_abs_vec(const Matrix& v)
{
    Matrix ret;
    details::eval_max_abs_vec::make<const Matrix&>(v,ret);
    return ret;
};

mat_int2 matcl::min2_vec(const Matrix& v)
{
    mat_int2 ret;
    details::eval_min2_vec::make<const Matrix&>(v,ret);
    return ret;
};

mat_int2 matcl::max2_vec(const Matrix& v)
{
    mat_int2 ret;
    details::eval_max2_vec::make<const Matrix&>(v,ret);
    return ret;
};

mat_int2 matcl::min_abs2_vec(const Matrix& v)
{
    mat_int2 ret;
    details::eval_min_abs2_vec::make<const Matrix&>(v,ret);
    return ret;
};

mat_int2 matcl::max_abs2_vec(const Matrix& v)
{
    mat_int2 ret;
    details::eval_max_abs2_vec::make<const Matrix&>(v,ret);
    return ret;
};

};
