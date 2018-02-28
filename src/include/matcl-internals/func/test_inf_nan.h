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

#include "matcl-internals/func/test_functor.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace details
{

struct test_inf_nan
{
    static const bool test_complex = true;

    template<class T>
    bool eval(const T& val)	
    { 
        return raw::details::isfinite_helper<T>::eval(val); 
    };

    bool eval(const Object& val)	
    { 
        return (bool)raw::details::isfinite_helper<Object>::eval(val); 
    };
};

struct test_nan
{
    static const bool test_complex = true;

    template<class T>
    bool eval(const T& val)	
    { 
        return raw::details::isnan_helper<T>::eval(val) == false; 
    };
    
    bool eval(const Object& val)	
    { 
        return (bool)raw::details::isnan_helper<Object>::eval(val) == false; 
    };
};

struct test_inf
{
    static const bool test_complex = true;

    template<class T>
    bool eval(const T& val)	
    { 
        return raw::details::isinf_helper<T>::eval(val) == false; 
    };
    
    bool eval(const Object& val)	
    { 
        return (bool)raw::details::isinf_helper<Object>::eval(val) == false; 
    };
};

// Test if value is not a nan and remember if it encountered an Inf
struct test_nan_remember_inf
{
    public:
        test_nan_remember_inf() :m_has_no_inf(1){};
        static const bool test_complex = true;

        template<class T>
        bool eval(const T& val)	
        {
            if (raw::details::isfinite_helper<T>::eval(val) == true)
            {
                return true;
            }
            else
            {
                m_has_no_inf = m_has_no_inf && (raw::details::isinf_helper<T>::eval(val) == false);
                return raw::details::isnan_helper<T>::eval(val) == false; 
            };
        };

        bool eval(const Object& val)	
        { 
            if ((bool)raw::details::isfinite_helper<Object>::eval(val) == true)
            {
                return true;
            }
            else
            {
                m_has_no_inf = m_has_no_inf && ((bool)raw::details::isinf_helper<Object>::eval(val) == false);
                return (bool)raw::details::isnan_helper<Object>::eval(val) == false; 
            }
        };
        
        // Get the recent has_no_inf; return false if Inf was encountered
        // since last visit
        bool get_has_no_inf() 
        {
            return m_has_no_inf;
        }

    private:
        bool m_has_no_inf;
};

template<class test_func>
struct test_matrix_visitor : public extract_type_switch<bool,test_matrix_visitor<test_func>,true>
{
    template<class T>
    static bool eval(const Matrix&, const T& mat, test_func& func)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return test_range<V, S, test_func>::eval(mat, func);
    };

    template<class T>
    static bool eval_scalar(const Matrix& handle, const T& mat, test_func& func)
    {
        (void)handle;
        return func.eval(mat);
    };
};

struct test_matrix
{
    template<class test_func>
    static bool eval(const matcl::Matrix& mat, test_func& f)
    {
        return test_matrix_visitor<test_func>::template make<const Matrix&>(mat,f);
    };
};

};}