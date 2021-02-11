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

#pragma once

#include "matcl-dynamic/matcl_dynamic.h"
#include "scalar.h"
#include <set>

namespace matcl { namespace test
{

namespace mdy = matcl::dynamic;

class test_dynamic
{
    private:
        using type2     = std::pair<mdy::Type, mdy::Type>;
        using type2_set = std::set<type2>;

    private:
        bool    test_predefined_unify_impl();
        bool    test_predefined_func_impl(const dynamic::function_name& func, bool promote_int);
        bool    test_object_impl();
        bool    test_one_impl();
        template<class T>
        double  test_one_type(bool has_one);
        bool    test_object_type_impl();
        bool    test_special_types_impl();
        bool    test_cons_impl();
        double  test_cons_1(mdy::Type t1, mdy::Type t2, const type2_set& missing, 
                        const type2_set& missing_as);
        bool    test_cons_val_impl();
        double  test_cons_val_1(const Scalar& s1, const Scalar& s2, Integer code);
        bool    test_object_func_impl();
        double  test_func_1(mdy::Type t1);
        bool    test_object_func_val_impl();
        double  test_func_1_val(const Scalar& s1, Integer code);
        bool    test_func_compare_impl();
        double  test_compare_1(mdy::Type t1, mdy::Type t2, const type2_set& missing);
        bool    test_func_uminus_impl();
        double  test_func_uminus_1(mdy::Type t1);
        bool    test_func_reim_impl();
        double  test_func_reim_1(mdy::Type t1);

        void    test_basic();
        void    test_object();
        void    test_one();
        void    test_object_type();
        void    test_special_types();
        void    test_cons();
        void    test_cons_val();
        void    test_object_func();
        void    test_object_func_val();
        void    test_func_compare();
        void    test_func_uminus();
        void    test_func_reim();
        void    test_predefined_unify();
        void    test_predefined_func(const dynamic::function_name& func, bool promote_int);
        void    test_compite_func();

        template<class T>
        void    test_compile_func_unary();
        template<class T1, class T2>
        void    test_compile_func_binary();
        template<class T, class T2>
        void    test_compile_func_binary_notcompl();
        template<class T, class T2>
        void    test_compile_func_binary_notmp();
        template<class T, class T2>
        void    test_compile_func_binary_str();
        template<class T>
        void    test_compile_func_unary_str();
        template<class T>
        void    test_compile_func_unary_notcompl();
        template<class T>
        void    test_compile_func_unary_notmp();

        std::set<type2> get_missing_cons();
        std::set<type2> get_missing_assign();

    public:
        void    make();
};

}};
