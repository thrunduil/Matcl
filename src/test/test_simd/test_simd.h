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

#include "test_simd_config.h"
#include "matcl-core/config.h"
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-simd/simd.h"

#include <string>

namespace matcl { namespace test
{

void test_performance_real();
void test_values_real();

class test_simd
{
    private:
        std::string m_instr_tag;
        bool        m_test_values;

    public:
        test_simd(bool test_values);

        void    make_ternary();
        void    make_binary();
        void    make_unary();

    private:
        int     get_N() const;
        int     get_M() const;

        template<class T>
        bool    test_equal(int size, const T* res, const T* res_gen, double max_dist, double& dist);

        template<class T>
        bool    test_equal(const T& res, const T& res_gen, double max_dist, double& dist);

        template<class T>
        void    test_functions();

        template<class T>
        void    test_functions_bin();

        template<class T>
        void    test_functions_3();

        template<class T, class Func>
        void    test_function(formatted_disp& fd, int size, int n_rep, const T* in, 
                    T* out, T* out_gen);

        template<class T, class Func>
        void    test_function_bin(formatted_disp& fd, int size, int n_rep, const T* in_1, 
                    const T* in_2, T* out, T* out_gen);

        template<class T, class Func>
        void    test_function_3(formatted_disp& fd, int size, int n_rep, const T* in_1, 
                    const T* in_2, const T* in_3, T* out, T* out_gen);

        template<class T, class Func>
        void    test_function_block(formatted_disp& fd, int size, int n_rep, const T* in, 
                    T* out, T* out_gen);

        template<class T, class Simd_type, class Func>
        double  test_function_simd(int size, int n_rep, const T* in, T* out);

        template<class T, class Simd_type, class Func>
        double  test_function_bin_simd(int size, int n_rep, const T* in1, const T* in2, T* out);

        template<class T, class Simd_type, class Func>
        double  test_function_3_simd(int size, int n_rep, const T* in1, const T* in2, const T* in3, T* out);

        template<class T, class Func>
        double  test_function_generic(int size, int n_rep, const T* in, T* out);
};

}}
