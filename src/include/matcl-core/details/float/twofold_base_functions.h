/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/config.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace details
{

template<class Float_type>
struct func_sqrt{};

template<class Float_type>
struct func_fms_a{};

template<class Float_type>
struct func_abs{};

//-----------------------------------------------------------------------
//                      func_signbit
//-----------------------------------------------------------------------
template<class T, int Bits, class Tag>
struct func_abs<matcl::simd::simd<T, Bits, Tag>>
{
    using simd_type = matcl::simd::simd<T, Bits, Tag>;

    force_inline
    static void eval(const simd_type& val, const simd_type& err,
                     simd_type& abs_val, simd_type& abs_err)
    {
        simd_type sign  = signbit_base(val);
        abs_val = bitwise_andnot(sign, val);
        abs_err = bitwise_andnot(sign, err);
    };
};

template<>
struct func_abs<float>
{
    force_inline
    static void eval(const float& val, const float& err,
                     float& abs_val, float& abs_err)
    {
        namespace mrds = matcl::raw::details::scal_func;

        // we cannot use comparison operator due to existence of -0.0
        bool sign       = mrds::signbit(val);

        if (sign == true)
        {
            abs_val = -val;
            abs_err = -err;
        }
        else
        {
            abs_val = val;
            abs_err = err;
        };
    };
};

template<>
struct func_abs<double>
{
    force_inline
    static void eval(const double& val, const double& err,
                     double& abs_val, double& abs_err)
    {
        namespace mrds = matcl::raw::details::scal_func;

        // we cannot use comparison operator due to existence of -0.0
        bool sign       = mrds::signbit(val);

        if (sign == true)
        {
            abs_val = -val;
            abs_err = -err;
        }
        else
        {
            abs_val = val;
            abs_err = err;
        };
    };
};

//-----------------------------------------------------------------------
//                      func_sqrt
//-----------------------------------------------------------------------
template<class T, int Bits, class Tag>
struct func_sqrt<matcl::simd::simd<T, Bits, Tag>>
{
    using simd_type = matcl::simd::simd<T, Bits, Tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd::sqrt(x);
    }
};

template<>
struct func_sqrt<float>
{
    force_inline
    static float eval(float x)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::sqrt(x);
    };
};

template<>
struct func_sqrt<double>
{
    force_inline
    static double eval(double x)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::sqrt(x);
    };
};

//-----------------------------------------------------------------------
//                      func_fms_a
//-----------------------------------------------------------------------
template<class T, int Bits, class Tag>
struct func_fms_a<matcl::simd::simd<T, Bits, Tag>>
{
    using simd_type = matcl::simd::simd<T, Bits, Tag>;

    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd::fms_a(x, y, z);
    };
};

template<>
struct func_fms_a<float>
{
    force_inline
    static float eval(const float& x, const float& y, const float& z)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::fms_a(x, y, z);
    };
};

template<>
struct func_fms_a<double>
{
    force_inline
    static double eval(const double& x, const double& y, const double& z)
    {
        namespace mrds = matcl::raw::details::scal_func;
        return mrds::fms_a(x, y, z);
    };
};

//-----------------------------------------------------------------------
//                      func_float_distance
//-----------------------------------------------------------------------
template<class Float_type, bool Is_simd>
struct func_float_distance
{};

template<class Float_type>
struct func_float_distance<Float_type, false>
{
    MATCL_CORE_EXPORT
    static Float_type
    eval(const twofold<Float_type>& x, const twofold<Float_type>& y);
};

template<class Float_type>
struct func_float_distance<Float_type, true>
{
    static Float_type
    eval(const twofold<Float_type>& x, const twofold<Float_type>& y)
    {
        static const int 
        vector_size         = Float_type::vector_size;

        using value_type    = typename Float_type::value_type;
        using twofold_base  = twofold<value_type>;

        Float_type res;

        const value_type* ptr_xv    = x.value.get_raw_ptr();
        const value_type* ptr_xe    = x.error.get_raw_ptr();
        const value_type* ptr_yv    = y.value.get_raw_ptr();
        const value_type* ptr_ye    = y.error.get_raw_ptr();

        value_type* ptr_res         = res.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
        {
            twofold_base xl = twofold_base(ptr_xv[i], ptr_xe[i]);
            twofold_base yl = twofold_base(ptr_yv[i], ptr_ye[i]);

            value_type dist = matcl::float_distance(xl, yl);

            ptr_res[i]      = dist;
        };

        return res;
    };
};

//-----------------------------------------------------------------------
//                      func_eps
//-----------------------------------------------------------------------
template<class Float_type, bool Is_simd>
struct func_eps
{};

template<>
struct func_eps<float, false>
{
    static float eval(const twofold<float>& x)
    {
        namespace mrds = matcl::raw::details::scal_func;

        float e = mrds::eps(x.value) * matcl::constants::eps<float>() * 0.5f;
        return e;
    };
};

template<>
struct func_eps<double, false>
{
    static double eval(const twofold<double>& x)
    {
        namespace mrds = matcl::raw::details::scal_func;

        double e = mrds::eps(x.value) * matcl::constants::eps() * 0.5;
        return e;
    };
};

template<class Float_type>
struct func_eps<Float_type, true>
{
    static Float_type eval(const twofold<Float_type>& x)
    {
        static const int 
        vector_size         = Float_type::vector_size;

        using value_type    = typename Float_type::value_type;
        using twofold_base  = twofold<value_type>;

        Float_type res;

        const value_type* ptr_xv    = x.value.get_raw_ptr();
        const value_type* ptr_xe    = x.error.get_raw_ptr();

        value_type* ptr_res         = res.get_raw_ptr();

        for (int i = 0; i < vector_size; ++i)
        {
            twofold_base xl = twofold_base(ptr_xv[i], ptr_xe[i]);
            value_type val  = matcl::eps(xl);

            ptr_res[i]      = val;
        };

        return res;
    };
};

//-----------------------------------------------------------------------
//                      func_save
//-----------------------------------------------------------------------
template<class Float_type, bool Is_simd>
struct func_save
{};

template<class Float_type>
struct func_save<Float_type, false>
{
    MATCL_CORE_EXPORT
    static void eval(std::ostream& os, const twofold<Float_type>& x);
};

template<class Float_type>
struct func_save<Float_type, true>
{
    static void eval(std::ostream& os, const twofold<Float_type>& x)
    {
        os << "{";
        os << x.value;
        os << ", ";
        os << x.error;
        os << "}";
    };
};

//-----------------------------------------------------------------------
//                      func_load
//-----------------------------------------------------------------------

template<class Float_type, bool Is_simd>
struct func_load
{};

template<class Float_type>
struct func_load<Float_type, false>
{
    MATCL_CORE_EXPORT
    static void eval(std::istream& is, twofold<Float_type>& x);
};

template<class Float_type>
struct func_load<Float_type, true>
{
    static void eval(std::istream& is, twofold<Float_type>& v)
    {
        using simd_type     = Float_type;
        using value_type    = typename simd_type::value_type;
  
        char c  = 0;

        // consume whitespaces
        while (is)
        {
            is.get(c);

            if (c != ' ' && c != '\t'  && c != '\n')
                break;
        }

        value_type nan  = std::numeric_limits<value_type>::quiet_NaN();
        simd_type val   = simd_type(nan);
        simd_type err   = simd_type(nan);
        v               = twofold<Float_type>(val, err);

        if (is.good() == false)
            return;

        if (c != '{')
        {
            // this is an error
            is.setstate(std::ios::failbit);
            return;
        }

        char sep    = ' ';
        char fin    = ' ';        
    
        is >> val;
        is >> sep;
        is >> err;
        is >> fin;

        if (sep != ',' || fin != '}')
        {
            // this is an error
            is.setstate(std::ios::failbit);
            return;
        };

        v = twofold<Float_type>(val, err);
        return;
    };
};

}}


