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

#include "matcl-matrep/func/raw/mvgen.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-internals/base/utils.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-internals/func/converter.h"

#pragma warning (disable:4127)

namespace matcl { namespace raw
{

namespace md = matcl::details;
namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;

integer_dense mr::irange(Integer s, Integer i, Integer e)
{
    if ((i == 0) || ((e - s) < 0 && i > 0) || ((e - s) > 0 && i < 0))
        return integer_dense(ti::ti_empty(),1,0);

    Integer j, n = (e - s) / i + 1, jj;

    integer_dense res(ti::ti_empty(),1,n);
    Integer * ptr = res.ptr();

    for (j = 0, jj = s; j < n; ++j, jj += i)
        *(ptr++) = jj;

    return res;
}

real_dense mr::range(Real s, Real i, Real e)
{
    if (i == 0 || (e-s)*i < 0)
        return real_dense(ti::ti_empty(),1,0);

    Real en = s + std::floor((e - s) / i) * i;
    Integer ii, n = icast(round(1. + (en - s) / i)); 

    real_dense res(ti::ti_empty(),1,n);

    Real * ptr = res.ptr();
    for (ii = 0; ii < n; ++ii)
        *(ptr++) = s + i * Real(ii);

    return res;
}

float_dense mr::frange(Float s, Float i, Float e)
{
    if (i == 0 || (e-s)*i < 0)
        return float_dense(ti::ti_empty(),1,0);

    Float en = s + std::floor((e - s) / i) * i;
    Integer ii, n = icast(round(1. + (en - s) / i)); 

    float_dense res(ti::ti_empty(),1,n);

    Float * ptr = res.ptr();
    for (ii = 0; ii < n; ++ii)
        *(ptr++) = s + i * Float(ii);

    return res;
}

real_dense mr::linspace(Real s, Real e, Integer n)
{
    if (n <= 0)
        return real_dense(ti::ti_empty(), 1, 0);

    if (n == 1)
        return real_dense(ti::ti_empty(),e,1, 1);

    Integer i;
    real_dense res(ti::ti_empty(),1,n);

    Real * ptr  = res.ptr();
    for (i = 1; i <= n; ++i)
        *(ptr++) = s + (e - s) * Real(i-1)/Real(n - 1);

    return res;
}

float_dense mr::flinspace(Float s, Float e, Integer n)
{
    if (n <= 0)
        return float_dense(ti::ti_empty(), 1, 0);

    if (n == 1)
        return float_dense(ti::ti_empty(),e,1, 1);

    Integer i;
    float_dense res(ti::ti_empty(),1,n);

    Float * ptr  = res.ptr();
    for (i = 1; i <= n; ++i)
        *(ptr++) = s + (e - s) * Float(i-1)/Float(n - 1);

    return res;
}

real_dense mr::logspace(Real s, Real e, Integer n)
{
    if (n <= 0)
        return raw::real_dense(ti::ti_empty(), 1, 0);

    if (n == 1)
    {
        //Matlab::logspace has strange, different semantics for e == pi.
        //return (e == constants::pi()) ? raw::real_dense(ti::ti_empty(),constants::pi(),1, 1) :
        //                     raw::real_dense(ti::ti_empty(),matcl::pow_nc(10., e),1, 1);
        
        Real val    = std::pow(10., e);
        return raw::real_dense(ti::ti_empty(),val,1, 1);
    };

    Integer i;
    real_dense res(ti::ti_empty(),1,n);

    //Matlab::logspace has strange, different semantics for e == pi.

    Real * ptr = res.ptr();
    for (i = 1; i <= n; ++i)
    {
        Real val    = s + (e - s) * Real(i-1)/Real(n - 1);
        val         = std::pow(10.0, val);
        *(ptr++)    = val;
    };

    return res;
}

float_dense mr::flogspace(Float s, Float e, Integer n)
{
    if (n <= 0)
        return raw::float_dense(ti::ti_empty(), 1, 0);

    if (n == 1)
    {
        Float val   = std::pow(10.f, e);
        return raw::float_dense(ti::ti_empty(),val,1, 1);
    };

    Integer i;
    float_dense res(ti::ti_empty(),1,n);

    Float * ptr     = res.ptr();
    for (i = 1; i <= n; ++i)
    {
        Float val   = s + (e - s) * Float(i-1)/Float(n - 1);
        val         = std::pow(10.0f, val);
        *(ptr++)    = val;
    };

    return res;
}

real_dense mr::eye(Integer r, Integer c)
{
    real_dense res(ti::ti_empty(),0.0, r, c);

    Integer i, n = (r < c) ? r : c, rp1 = r + 1;

    if ( n == 0)
        return res;

    Real * ptr = res.ptr();
    for (i = 0; i < n; ++i)
    {
        *ptr = 1.;
        ptr += rp1;
    };

    if (r == c)
        res.get_struct().set(predefined_struct_type::id);
    else
        res.get_struct().set(predefined_struct_type::diag);

    return res;
}

integer_dense mr::ieye(Integer r, Integer c)
{
    integer_dense res(ti::ti_empty(),Integer(), r, c);

    Integer i, n = (r < c) ? r : c, rp1 = r + 1;

    if ( n == 0)
        return res;

    Integer * ptr = res.ptr();

    for (i = 0; i < n; ++i)
    {
        *ptr = 1;
        ptr += rp1;
    };

    if (r == c)
        res.get_struct().set(predefined_struct_type::id);
    else
        res.get_struct().set(predefined_struct_type::diag);

    return res;
}

complex_dense mr::ceye(Integer r, Integer c)
{
    complex_dense res(ti::ti_empty(), 0.0, r, c);

    Integer n = (r < c) ? r : c, rp1 = r + 1;

    if ( n == 0)
        return res;

    Complex * ptr = res.ptr();

    for (Integer i = 0; i < n; ++i)
    {
        *ptr = 1.;
        ptr += rp1;
    };

    if (r == c)
        res.get_struct().set(predefined_struct_type::id);
    else
        res.get_struct().set(predefined_struct_type::diag);

    return res;
}

float_dense mr::feye(Integer r, Integer c)
{
    float_dense res(ti::ti_empty(), 0.0f, r, c);

    Integer n = (r < c) ? r : c, rp1 = r + 1;

    if ( n == 0)
        return res;

    Float * ptr = res.ptr();

    for (Integer i = 0; i < n; ++i)
    {
        *ptr = 1.f;
        ptr += rp1;
    };

    if (r == c)
        res.get_struct().set(predefined_struct_type::id);
    else
        res.get_struct().set(predefined_struct_type::diag);

    return res;
}

float_complex_dense mr::fceye(Integer r, Integer c)
{
    float_complex_dense res(ti::ti_empty(), 0.0f, r, c);

    Integer n = (r < c) ? r : c, rp1 = r + 1;

    if ( n == 0)
        return res;

    Float_complex * ptr = res.ptr();

    for (Integer i = 0; i < n; ++i)
    {
        *ptr = 1.f;
        ptr += rp1;
    };

    if (r == c)
        res.get_struct().set(predefined_struct_type::id);
    else
        res.get_struct().set(predefined_struct_type::diag);

    return res;
}

object_dense mr::eye(ti::ti_object ti,Integer r, Integer c)
{
    ti::ti_object ret_ti = ti;

    object_dense res(ret_ti,md::default_value<Object>(ti), r, c);

    Integer i, n = (r < c) ? r : c, rp1 = r + 1;

    if ( n == 0)
        return res;

    Object one = md::one_value<Object>(ret_ti);

    Object * ptr = res.ptr();
    for (i = 0; i < n; ++i)
    {
        mrd::reset_helper(*ptr, one);
        ptr += rp1;
    };

    if (r == c)
        res.get_struct().set(predefined_struct_type::id);
    else
        res.get_struct().set(predefined_struct_type::diag);

    return res;
}

template<class V>
struct speye_helper
{
    using SM = raw::Matrix<V,struct_sparse>;
    static SM eval(ti::ti_type<V> ti, Integer r, Integer c)
    {
        Integer m = (r < c) ? r : c;

        details::sparse_ccs<V> tmp(ti, r, c, m);
        if (m == 0)
            return sparse_matrix_base<V>(tmp);

        V one = md::one_value<V>(ti);

        V * ptr_x = tmp.ptr_x();
        for (Integer i = 0; i < m; ++i)
            mrd::reset_helper(*(ptr_x++), one);

        Integer * ptr_c = tmp.ptr_c();
        Integer * ptr_r = tmp.ptr_r();
        Integer i       = 0;
        for (; i < m; ++i)
        {
            *(ptr_c++) = i;
            *(ptr_r++) = i;
        }

        for (i = m; i <= c; ++i)
            *(ptr_c++) = m;

        SM out = sparse_matrix_base<V>(tmp);

        if (r == c)
            out.get_struct().set(predefined_struct_type::id);
        else
            out.get_struct().set(predefined_struct_type::diag);

        return out;
    };
};

real_sparse mr::speye(Integer r, Integer c)
{
    return speye_helper<Real>::eval(ti::ti_empty(),r,c);
}

real_sparse mr::speye(Integer n)
{
    return speye(n,n);
}

integer_sparse mr::ispeye(Integer r, Integer c)
{
    return speye_helper<Integer>::eval(ti::ti_empty(),r,c);
}

integer_sparse mr::ispeye(Integer r)
{
    return ispeye(r,r);
}

complex_sparse mr::cspeye(Integer r, Integer c)
{
    return speye_helper<Complex>::eval(ti::ti_empty(),r,c);
}

complex_sparse mr::cspeye(Integer r)
{
    return cspeye(r,r);
};

float_sparse mr::fspeye(Integer r, Integer c)
{
    return speye_helper<Float>::eval(ti::ti_empty(),r,c);
}

float_sparse mr::fspeye(Integer r)
{
    return fspeye(r,r);
};

float_complex_sparse mr::fcspeye(Integer r, Integer c)
{
    return speye_helper<Float_complex>::eval(ti::ti_empty(),r,c);
}

float_complex_sparse mr::fcspeye(Integer r)
{
    return fcspeye(r,r);
};

object_sparse mr::speye(ti::ti_object ti, Integer r, Integer c)
{
    return speye_helper<Object>::eval(ti,r,c);
}

object_sparse mr::speye(ti::ti_object ti, Integer n)
{
    return mr::speye(ti,n,n);
}

template<class V>
struct beye_helper
{
    using BM = raw::Matrix<V,struct_banded>;

    static BM eval(ti::ti_type<V> ti, Integer m, Integer n, Integer fd, Integer ld)
    {
        V Z = md::default_value<V>(ti);
        BM out(ti, Z, m, n, fd, ld);

        Integer out_size = std::min(m,n);
    
        if (out_size == 0)
            return out;

        V* ptr      = out.rep_ptr() + out.first_elem_diag(0);
        Integer ldo = out.ld();
        V one       = md::one_value<V>(ti);

        for (Integer i = 0; i < out_size; ++i, ptr+= ldo)
            mrd::reset_helper(*ptr,one);

        if (m == n)
            out.get_struct().set(predefined_struct_type::id);
        else
            out.get_struct().set(predefined_struct_type::diag);

        return out;
    };
};

real_band mr::beye(Integer m, Integer n, Integer fd, Integer ld)
{
    return beye_helper<Real>::eval(ti::ti_empty(),m,n,fd,ld);
};

real_band mr::beye(Integer n, Integer fd, Integer ld)
{
    return mr::beye(n,n,fd,ld);
};

object_band mr::beye(ti::ti_object ti,Integer m, Integer n, Integer fd, Integer ld)
{
    return beye_helper<Object>::eval(ti,m,n,fd,ld);
};

object_band mr::beye(ti::ti_object ti,Integer n, Integer fd, Integer ld)
{
    return mr::beye(ti,n,n,fd,ld);
};

integer_band mr::ibeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return beye_helper<Integer>::eval(ti::ti_empty(),m,n,fd,ld);
};

integer_band mr::ibeye(Integer n, Integer fd, Integer ld)
{
    return mr::ibeye(n,n,fd,ld);
};

complex_band mr::cbeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return beye_helper<Complex>::eval(ti::ti_empty(),m,n,fd,ld);
};

complex_band mr::cbeye(Integer n, Integer fd, Integer ld)
{
    return mr::cbeye(n,n,fd,ld);
};

float_band mr::fbeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return beye_helper<Float>::eval(ti::ti_empty(),m,n,fd,ld);
};

float_band mr::fbeye(Integer n, Integer fd, Integer ld)
{
    return mr::fbeye(n,n,fd,ld);
};

float_complex_band mr::fcbeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return beye_helper<Float_complex>::eval(ti::ti_empty(),m,n,fd,ld);
};

float_complex_band mr::fcbeye(Integer n, Integer fd, Integer ld)
{
    return mr::fcbeye(n,n,fd,ld);
};

template<class ret_type, class in_type>
struct diag_helper
{
    using value_type_in     = typename in_type::value_type;
    using value_type_ret    = typename ret_type::value_type;

    static ret_type eval(ti::ti_type<value_type_ret> ret_ti, const in_type& v0, Integer d)
    { 
        in_type v = v0.make_explicit();

        value_type_ret Z = md::default_value<value_type_ret>(ret_ti);

        Integer n = v.size();
        error::check_diag_arg(v.rows(), v.cols());

        if (n == 0) 
            return ret_type(ret_ti);

        Integer r, c, st;
        if (d >= 0)
        {
            r   = n + d; 
            c   = n + d;
            st  = imult(d,r);
        }
        else
        {
            r   = n - d;
            c   = n - d;
            st  = - d;
        };
 
        ret_type res(ret_ti, Z, r, c);

        if (n==0)
            return res;

        value_type_ret* ptr = res.ptr() + st;
        const value_type_in* ptr_v = v.ptr();
        Integer res_ld  = res.ld();
        Integer ldv     = (v.rows() == 1)? v.ld() : 1;

        for (Integer i = 0; i < n; ++i)
        {
            mrd::reset_helper(*ptr,ptr_v[0]);
            ptr     += res_ld + 1;
            ptr_v   += ldv;
        };
 
        if (d == 0)     res.get_struct().set(predefined_struct_type::diag);
        else if (d > 0) res.get_struct().set(predefined_struct_type::triu);
        else            res.get_struct().set(predefined_struct_type::tril);

        return res;
    };
};

integer_dense mr::diag(const integer_dense &v, Integer d)
{
    return diag_helper<integer_dense,integer_dense>::eval(ti::ti_empty(),v,d);
};

real_dense mr::diag(const real_dense &v, Integer d)
{
    return diag_helper<real_dense,real_dense>::eval(ti::ti_empty(),v,d);
}

complex_dense mr::diag(const complex_dense &v, Integer d)
{
    return diag_helper<complex_dense,complex_dense>::eval(ti::ti_empty(),v,d);
}

float_dense mr::diag(const float_dense &v, Integer d)
{
    return diag_helper<float_dense,float_dense>::eval(ti::ti_empty(),v,d);
}

float_complex_dense mr::diag(const float_complex_dense &v, Integer d)
{
    return diag_helper<float_complex_dense,float_complex_dense>::eval(ti::ti_empty(),v,d);
}

object_dense mr::diag(const object_dense &v, Integer d)
{
    ti::ti_object ret_ti = v.get_type();
    return diag_helper<object_dense,object_dense>::eval(ret_ti,v,d);
};

template<class V>
struct bdiag_helper
{
    using DM = raw::Matrix<V,struct_dense>;
    using BM = raw::Matrix<V,struct_banded>;

    static BM eval(const DM& v0, Integer d)
    {
        DM v = v0.make_explicit();
        ti::ti_type<V> ret_ti = v.get_type();        

        error::check_bspdiag_1starg(v.rows(), v.cols());

        Integer dl = v.size();

        if (dl == 0) 
            return BM(ret_ti,0,0,0,0);

        Integer n   = (d >= 0) ? dl + d : dl - d;

        V Z = md::default_value<V>(ret_ti);
        BM res(ret_ti, Z, n, n, d, d);

        if (n == 0)
            return res;
        
        Integer ldr = res.ld();
        Integer ldv = (v.rows() == 1)? v.ld() : 1;

        V * ptr     = res.rep_ptr() + res.first_elem_diag(d);
        const V * ptr_in = v.ptr();

        for (Integer i = 0; i < dl; ++i)
        {
            mrd::reset_helper(*ptr,*ptr_in);
            ptr     += ldr;
            ptr_in  += ldv;
        };

        if (d == 0)     res.get_struct().set(predefined_struct_type::diag);
        else if (d > 0) res.get_struct().set(predefined_struct_type::triu);
        else            res.get_struct().set(predefined_struct_type::tril);

        return res;
    };
};

integer_band mr::bdiag(const integer_dense &v, Integer d)
{
    return bdiag_helper<Integer>::eval(v,d);
}

object_band mr::bdiag(const object_dense &v, Integer d)
{
    return bdiag_helper<Object>::eval(v,d);
}

integer_band mr::bdiag(const integer_sparse &v, Integer d)
{
    return bdiag(full(v),d);
};

real_band mr::bdiag(const real_dense &v, Integer d)
{
    return bdiag_helper<Real>::eval(v,d);
}

real_band mr::bdiag(const real_sparse &v, Integer d)
{
    return bdiag(full(v),d);
};

complex_band mr::bdiag(const complex_dense &v, Integer d)
{
    return bdiag_helper<Complex>::eval(v,d);
}

complex_band mr::bdiag(const complex_sparse &v, Integer d)
{
    return bdiag(full(v),d);
};

float_band mr::bdiag(const float_dense &v, Integer d)
{
    return bdiag_helper<Float>::eval(v,d);
}

float_band mr::bdiag(const float_sparse &v, Integer d)
{
    return bdiag(full(v),d);
};

float_complex_band mr::bdiag(const float_complex_dense &v, Integer d)
{
    return bdiag_helper<Float_complex>::eval(v,d);
}

float_complex_band mr::bdiag(const float_complex_sparse &v, Integer d)
{
    return bdiag(full(v),d);
};

template<class V>
struct spdiag_helper
{
    using DM = raw::Matrix<V,struct_dense>;
    using SM = raw::Matrix<V,struct_sparse>;

    static SM eval(const DM &v0, Integer d)
    {
        DM v = v0.make_explicit();

        error::check_bspdiag_1starg(v.rows(), v.cols());

        ti::ti_type<V> ret_ti = ti::get_ti(v);

        Integer dl = v.size();
        if (dl == 0)
        {
            details::sparse_ccs<V> dat(ret_ti,0, 0, 0);
            SM out = sparse_matrix_base<V>(dat);
            return out;
        };

        Integer ldv = (v.rows() == 1)? v.ld() : 1;
        Integer n;

        if (d >= 0)
            n = dl + d;
        else
            n = dl - d;

        details::sparse_ccs<V> dat(ret_ti,n, n, dl);

        if (d >= 0)
        {
            Integer * ptr_c = dat.ptr_c() + d+1;
            Integer * ptr_r = dat.ptr_r();
            V * ptr_x       = dat.ptr_x();
            const V * ptr_v = v.ptr();
            
            for (Integer j = d, i = 0; j < n; ++j, ++i)
            {
                *(ptr_c++) = i + 1;
                *(ptr_r++) = i;
                mrd::reset_helper(*(ptr_x++),*ptr_v);
                ptr_v += ldv;
            }
        }
        else
        {
            Integer * ptr_c = dat.ptr_c() + 1;
            Integer * ptr_r = dat.ptr_r();
            V * ptr_x       = dat.ptr_x();
            const V * ptr_v = v.ptr();

            Integer j, i;
            for (j = 1, i = -d; j <= dl; ++j, ++i)
            {
                *(ptr_c++) = j;
                *(ptr_r++) = i;
                mrd::reset_helper(*(ptr_x++),*ptr_v);
                ptr_v += ldv;
            }

            for (; j <= n; ++j, ++i)
                *(ptr_c++) = dl;
        }

        SM out = sparse_matrix_base<V>(dat);
        if (d == 0)     out.get_struct().set(predefined_struct_type::diag);
        else if (d > 0) out.get_struct().set(predefined_struct_type::triu);
        else            out.get_struct().set(predefined_struct_type::tril);

        return out;
    }
};

integer_sparse mr::spdiag(const integer_dense &v, Integer d)
{
    return spdiag_helper<Integer>::eval(v,d);
}

real_sparse mr::spdiag(const real_dense &v, Integer d)
{
    return spdiag_helper<Real>::eval(v,d);
}

complex_sparse mr::spdiag(const complex_dense &v, Integer d)
{
    return spdiag_helper<Complex>::eval(v,d);
}

float_sparse mr::spdiag(const float_dense &v, Integer d)
{
    return spdiag_helper<Float>::eval(v,d);
}

float_complex_sparse mr::spdiag(const float_complex_dense &v, Integer d)
{
    return spdiag_helper<Float_complex>::eval(v,d);
}

object_sparse mr::spdiag(const object_dense &v, Integer d)
{
    return spdiag_helper<Object>::eval(v,d);
}

template<class V>
struct bdiags_helper
{
    using DM = raw::Matrix<V,struct_dense>;
    using BM = raw::Matrix<V,struct_banded>;

    static BM eval(const DM& B, const integer_dense &d0, Integer m, Integer n)
    {
        integer_dense d = d0.make_explicit();

        error::check_bspdiags_2ndarg(d.rows(), d.cols());

        ti::ti_type<V> ret_ti = B.get_type();

        Integer mi = constants::max_int(), mx = constants::min_int();
        Integer s = d.size(), c = B.cols(), r = B.rows();

        if (s != c)
            throw error::bspdiags_nonconf();
        
        if (std::min(m,n) != r)
            throw error::bspdiags_nonconf();

        const Integer * ptr_d = d.ptr();
        Integer ldd = (d.rows() == 1)? d.ld() : 1;

        for (Integer k = 0; k < s; ++k)
        {
            if (*ptr_d < mi) mi = *ptr_d;
            if (*ptr_d > mx) mx = *ptr_d;

            ptr_d += ldd;
        }

        Integer fd = mi;
        Integer ld = mx;

        if (m < 0 || n < 0 || (n > 0 && ld + 1 > n) || (m > 0 && -fd + 1 > m) 
                || (n == 0 && ld > 0) || (m == 0 && fd < 0))
        {
            throw error::bspdiags_nonconf();
        };

        V Z = md::default_value<V>(ret_ti);

        BM res(ret_ti, Z, m, n, fd, ld);

        if (m == 0 || n == 0 )
            return res;

        V * ptr = res.rep_ptr();

        Integer ldr = res.ld();
        ptr_d       = d.ptr();

        for (Integer k = 0; k < s; ++k)
        {
            Integer diag    = ptr_d[0];
            ptr_d           += ldd;

            Integer fp      = res.first_elem_diag(diag);
            Integer nelem   = res.diag_length(diag);

            ptr             = res.rep_ptr() + fp;
            const V * ptr_B = B.ptr() + imult(k,B.ld());

            for (Integer i = 0; i < nelem; ++i)
            {
                mrd::reset_helper(*ptr,*(ptr_B++));
                ptr += ldr;
            };
        }

        return res;
    };
};

integer_band mr::bdiags(const integer_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    return bdiags_helper<Integer>::eval(B,d,m,n);
}

real_band mr::bdiags(const real_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    return bdiags_helper<Real>::eval(B,d,m,n);
}

complex_band mr::bdiags(const complex_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    return bdiags_helper<Complex>::eval(B,d,m,n);
}

float_band mr::bdiags(const float_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    return bdiags_helper<Float>::eval(B,d,m,n);
}

float_complex_band mr::bdiags(const float_complex_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    return bdiags_helper<Float_complex>::eval(B,d,m,n);
}

object_band mr::bdiags(const object_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    return bdiags_helper<Object>::eval(B,d,m,n);
}

integer_dense mr::diags(const integer_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    integer_band tmp = bdiags(B, d, m, n);
    return converter<integer_dense,integer_band>::eval(tmp);
}

real_dense mr::diags(const real_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    real_band tmp = bdiags(B, d, m, n);
    return converter<real_dense,real_band>::eval(tmp);
}

complex_dense mr::diags(const complex_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    complex_band tmp = bdiags(B, d, m, n);
    return converter<complex_dense,complex_band>::eval(tmp);
}

float_dense mr::diags(const float_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    float_band tmp = bdiags(B, d, m, n);
    return converter<float_dense,float_band>::eval(tmp);
}

float_complex_dense mr::diags(const float_complex_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    float_complex_band tmp = bdiags(B, d, m, n);
    return converter<float_complex_dense,float_complex_band>::eval(tmp);
}

object_dense mr::diags(const object_dense &B, const integer_dense &d,
              Integer m, Integer n)
{
    object_band tmp = bdiags(B, d, m, n);
    return converter<object_dense,object_band>::eval(tmp);
}

integer_sparse mr::spdiags(const integer_dense &B, const integer_dense &d, Integer m,
               Integer n)
{
    return sparse(bdiags(B, d, m, n));
}

real_sparse mr::spdiags(const real_dense &B, const integer_dense &d, Integer m,
               Integer n)
{
    return sparse(bdiags(B, d, m, n));
}

complex_sparse mr::spdiags(const complex_dense &B, const integer_dense &d, Integer m,
               Integer n)
{
    return sparse(bdiags(B, d, m, n));
}

float_sparse mr::spdiags(const float_dense &B, const integer_dense &d, Integer m,
               Integer n)
{
    return sparse(bdiags(B, d, m, n));
}

float_complex_sparse mr::spdiags(const float_complex_dense &B, const integer_dense &d, Integer m,
               Integer n)
{
    return sparse(bdiags(B, d, m, n));
}

object_sparse mr::spdiags(const object_dense &B, const integer_dense &d, Integer m,
               Integer n)
{
    return sparse(bdiags(B, d, m, n));
}

template<class Value, class Derived>
struct rand_helper
{
    using DM = mr::Matrix<Value, struct_dense>;

    static DM eval(const matcl::rand_state& rand_ptr, Integer r, Integer c)
    {
        DM res(ti::ti_empty(),r, c);

        Integer n = res.size();

        if (n == 0)
            return res;

        Value * ptr = res.ptr();
        for (Integer i = 0; i < n; ++i)
            ptr[i]  = Derived::rand(rand_ptr);

        return res;
    }
};

real_dense mr::rand(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Real, rand_impl>
    {
        static Real rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::rand(rand_ptr);
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
}

integer_dense mr::irand(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Integer, rand_impl>
    {
        static Integer rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::irand(rand_ptr);
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
}

complex_dense mr::crand(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Complex, rand_impl>
    {
        static Complex rand(const matcl::rand_state& rand_ptr)
        {
            return Complex(matcl::rand(rand_ptr), matcl::rand(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
};

float_dense mr::frand(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Float, rand_impl>
    {
        static Float rand(const matcl::rand_state& rand_ptr)
        {
            return static_cast<Float>(matcl::rand(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
};

float_complex_dense mr::fcrand(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Float_complex, rand_impl>
    {
        static Float_complex rand(const matcl::rand_state& rand_ptr)
        {
            return Float_complex(static_cast<Float>(matcl::rand(rand_ptr)),
                                 static_cast<Float>(matcl::rand(rand_ptr)));
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
};

real_dense mr::randn(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Real, rand_impl>
    {
        static Real rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::randn(rand_ptr);
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
}

complex_dense mr::crandn(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Complex, rand_impl>
    {
        static Complex rand(const matcl::rand_state& rand_ptr)
        {
            return Complex(matcl::randn(rand_ptr),matcl::randn(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
}

float_dense mr::frandn(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Float, rand_impl>
    {
        static Float rand(const matcl::rand_state& rand_ptr)
        {
            return static_cast<Float>(matcl::randn(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
}

float_complex_dense mr::fcrandn(const matcl::rand_state& rand_ptr, Integer r, Integer c)
{
    struct rand_impl : rand_helper<Float_complex, rand_impl>
    {
        static Float_complex rand(const matcl::rand_state& rand_ptr)
        {
            return Float_complex(static_cast<Float>(matcl::randn(rand_ptr)),
                                 static_cast<Float>(matcl::randn(rand_ptr)));
        };
    };

    return rand_impl::eval(rand_ptr,r,c);
}

integer_dense mr::randperm(const matcl::rand_state& rand_ptr, Integer N)
{
    error::check_randperm_arg(N);

    if (N == 0)
        return integer_dense(ti::ti_empty(),1,0);

    integer_dense res = irange(1, N);
    Integer* ptr_res = res.ptr();

    for (Integer i = 0, k = N; i < N; ++i, --k)
    {
        Integer sw  = ptr_res[i];
        Integer j   = i + icast(k * matcl::rand(rand_ptr));
        ptr_res[i]  = ptr_res[j];
        ptr_res[j]  = sw;
    }

    return res;
}

template<class Value, class Derived>
struct sprand_helper
{
    using SM = mr::Matrix<Value,struct_sparse>;
    using DM = mr::Matrix<Value,struct_dense>;

    static SM eval(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
    {
        if (d < 0.) d = 0.;
        if (d > 1.) d = 1.;

        Integer i, nnz = icast_c(d * Real(r) * Real(c));

        integer_dense ri(ti::ti_empty(),nnz,1), ci(ti::ti_empty(),nnz,1);
        Integer* ptr_ri = ri.ptr();
        Integer* ptr_ci = ci.ptr();

        DM xv(ti::ti_empty(),nnz,1);

        for (i = 0; i < nnz; ++i)
        {
            ptr_ri[i] = 1 + icast(matcl::rand(rand_ptr) * Real(r));
            ptr_ci[i] = 1 + icast(matcl::rand(rand_ptr) * Real(c));
        }

        SM tmp(ti::ti_empty(),ri.ptr(), ci.ptr(), xv.ptr(), r, c, nnz, nnz);

        Integer nz      = tmp.nnz();
        Value * ptr_x   = tmp.rep().ptr_x();

        for (i = 0; i < nz; ++i)
            ptr_x[i]    = Derived::rand(rand_ptr);

        return tmp;
    }
};

real_sparse mr::sprand(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Real,rand_impl>
    {
        static Real rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::rand(rand_ptr);
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

float_sparse mr::fsprand(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Float,rand_impl>
    {
        static Float rand(const matcl::rand_state& rand_ptr)
        {
            return static_cast<Float>(matcl::rand(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

complex_sparse mr::csprand(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Complex,rand_impl>
    {
        static Complex rand(const matcl::rand_state& rand_ptr)
        {
            return Complex(matcl::rand(rand_ptr), matcl::rand(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

float_complex_sparse mr::fcsprand(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Float_complex,rand_impl>
    {
        static Float_complex rand(const matcl::rand_state& rand_ptr)
        {
            return Float_complex(static_cast<Float>(matcl::rand(rand_ptr)),
                                 static_cast<Float>(matcl::rand(rand_ptr)));
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

integer_sparse mr::isprand(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Integer,rand_impl>
    {
        static Integer rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::irand(rand_ptr);
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

real_sparse mr::sprandn(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Real,rand_impl>
    {
        static Real rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::randn(rand_ptr);
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

float_sparse mr::fsprandn(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Float,rand_impl>
    {
        static Float rand(const matcl::rand_state& rand_ptr)
        {
            return static_cast<Float>(matcl::randn(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

complex_sparse mr::csprandn(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Complex,rand_impl>
    {
        static Complex rand(const matcl::rand_state& rand_ptr)
        {
            return Complex(matcl::randn(rand_ptr),matcl::randn(rand_ptr));
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

float_complex_sparse mr::fcsprandn(const matcl::rand_state& rand_ptr, Integer r, Integer c, Real d)
{
    struct rand_impl : sprand_helper<Float_complex,rand_impl>
    {
        static Float_complex rand(const matcl::rand_state& rand_ptr)
        {
            return Float_complex(static_cast<Float>(matcl::randn(rand_ptr)),
                                 static_cast<Float>(matcl::randn(rand_ptr)));
        };
    };

    return rand_impl::eval(rand_ptr, r, c, d);
}

};};