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
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/details/isa.h"
#include "matcl-core/details/type_codes.h"

namespace matcl
{

// manipulation of struct_type, value_type and mat_code
struct MATCL_CORE_EXPORT matrix_traits
{
    private:
        
        template<class mat_type,bool is_scal>
        struct mat_type_info_type_impl
        {
            using matrix_type   = mat_type;
            using value_type    = typename matrix_type::value_type;
            using struct_type   = typename matrix_type::struct_type;

            static const matcl::mat_code
            matrix_code         = details::type_to_code<mat_type>::value;

            static const matcl::value_code
            value_code          = details::value_to_code<value_type>::value;

            static const matcl::struct_code
            struct_code         = details::struct_to_code<struct_type>::value;
        };

        template<class mat_type>
        struct mat_type_info_type_impl<mat_type,true>
        {
            using matrix_type   = mat_type;
            using value_type    = mat_type;
            using struct_type   = struct_scalar;

            static const matcl::mat_code 
            matrix_code         = details::type_to_code<mat_type>::value;

            static const matcl::value_code
            value_code          = details::value_to_code<value_type>::value;

            static const matcl::struct_code
            struct_code         = details::struct_to_code<struct_type>::value;
        };


    // internal use only
    public:    
        template<class mat_type>
        struct mat_type_info_type : 
            public mat_type_info_type_impl<mat_type,details::is_scalar<mat_type>::value>
        {};

        template<class vt, class st>
        struct mat_type_info_type_2 : public mat_type_info_type<raw::Matrix<vt,st>>
        {};

    public:
        // get value code from matrix code in runtime
        static matcl::value_code    get_value_type(matcl::mat_code mt);

        // get structure code from matrix code in runtime
        static matcl::struct_code   get_struct_type(matcl::mat_code mt);
        
        // get matrix code from value code and structure code in runtime
        static matcl::mat_code      get_matrix_type(matcl::value_code vt, 
                                        matcl::struct_code st);

        // get value code and structure code from matrix code at compile time
        template<matcl::mat_code mt>
        struct mat_type_info_code
        {
            // internal use
            using matrix_type   = typename details::code_to_type<mt>::type;
        
            // internal use
            using struct_type   = typename matrix_type::struct_type;

            // type of stored elements (Integer, Float, Real, Complex, 
            // Float_complex, or Object)
            using value_type    = typename matrix_type::value_type;

            // matrix code 
            static const matcl::mat_code
            matrix_code         = mt;
            
            // value code
            static const matcl::value_code
            value_code          = details::value_to_code<value_type>::value;
            
            // structure code
            static const matcl::struct_code
            struct_code         = details::struct_to_code<struct_type>::value;
        };

        // get matrix code from value code and struct code at compile time
        template<matcl::value_code vt, matcl::struct_code st>
        struct mat_type_info_code_2
        {
            // type of stored elements (Integer, Float, Real, Complex, 
            // Float_complex, or Object)
            using value_type    = typename details::code_to_value<vt>::type;
            
            // internal use
            using struct_type   = typename details::code_to_struct<st>::type;
            
            // internal use
            using matrix_type   = raw::Matrix<value_type,struct_type>;

            // matrix code 
            static const matcl::mat_code
            matrix_code         = details::type_to_code<matrix_type>::value;
            
            // value code
            static const matcl::value_code
            value_code          = vt;
            
            // structure code
            static const matcl::struct_code 
            struct_code         = st;
        };

        // get value code from value type; V must be valid scalar type
        // i.e. Integer, Float, Real, Complex, Float_complex, or Object
        template<class V>
        struct value_code
        {
            // value code of given value type
            static const matcl::value_code
            value               = details::value_to_code<V>::value;
        };

        // get value code of the type T, such that conversions form T1 to T and from
        // T2 to T are best possible, where v1, v2 are value codes of T1 and T2
        static matcl::value_code unify_value_types(matcl::value_code v1, 
                                    matcl::value_code v2);

        // get value code of the real type of the type T, where v1 is value codes of T
        static matcl::value_code real_value_type(matcl::value_code v1);

        // get value code of the complex type of the type T, where v1 is value codes of T
        static matcl::value_code complex_value_type(matcl::value_code v1);

        // return true if value code v represents real floating point value or Integer
        static bool             is_real(matcl::value_code v);

        // return true if value code v represents floating point value (real or complex)
        static bool             is_float(matcl::value_code v);

        // return true if value code v represents real floating point value
        static bool             is_float_real(matcl::value_code v1);

        // return true if value code v represents complex floating point value
        static bool             is_float_complex(matcl::value_code v1);

        // return true if value code v represents Float or Float_complex
        static bool             is_single_precision(matcl::value_code vc);

        // convert given value code to equivalent value code representing
        // double precision scalars; Integer is converted to Real, Object is
        // not changed; if v represents a complex scalar, then returned value
        // code also represents a complex scalar
        static matcl::value_code double_precision(matcl::value_code v);

        // convert given value code to equivalent value code representing
        // single precision scalars; Integer is converted to Float, Object is
        // not changed; if v represents a complex scalar, then returned value
        // code also represents a complex scalar
        static matcl::value_code single_precision(matcl::value_code v);
};

};