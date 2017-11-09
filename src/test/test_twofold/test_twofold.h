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
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-core/float/twofold.h"

#include <string>

namespace matcl { namespace test
{

void test_fma();
void test_double();
void test_io();
void test_error();

class test_twofold
{
    private:
        bool    m_normalized;

    public:
        test_twofold(bool normalized);

        void    make_fma();
        void    make_binary();
        void    make_unary();
        void    make_io();

        void    make_error();

    private:
        int     get_N() const;
        int     get_M() const;
        
        bool    make_error(const twofold& x, int mult);

        template<class T>
        std::string get_header() const;

        template<class T>
        void    test_functions();

        template<class T>
        void    test_functions_io();

        template<class T>
        void    test_functions_bin();

        template<class T>
        void    test_functions_fma();

        template<class T, class T2, class TMP>
        void    test_functions_bin(formatted_disp& fd, int size, int n_rep, T2* in_1, 
                    T2* in_2, T2* out, TMP* out_gen);

        template<class T, class T2, class TMP>
        void    test_functions_sum_bin(formatted_disp& fd, int size, int n_rep, T2* in_1, 
                    T2* in_2, T2* out, TMP* out_gen);

        template<class T, class T2, class TMP, class Func>
        void    test_function(formatted_disp& fd, int size, int n_rep, const T2* in, 
                    T2* out, TMP* out_gen, double max_dist);

        template<class T, class T2, class Func>
        void    test_function_io(formatted_disp& fd, int size, int n_rep, const T2* in, 
                    T2* out, double max_dist);

        template<class T, class T2, class TMP, class Func>
        void    test_function_bin(formatted_disp& fd, int size, int n_rep, const T2* in_1, 
                    const T2* in_2, T2* out, TMP* out_gen, double max_dist);

        template<class T, class TMP, class Func>
        void    test_function_fma(formatted_disp& fd, int size, int n_rep, const T* in_1, 
                    const T* in_2, const T* in_3, T* out, TMP* out_gen, double max_dist);

        template<class T2, class TMP>
        bool    test_equal(int size, const T2* res, const TMP* res_gen, double max_dist, double& dist);

        template<class T2>
        bool    test_equal(int size, const T2* res, const T2* res_gen, double max_dist, double& dist);

        template<class T2, class TMP>
        bool    test_equal(const T2& res, const TMP& res_gen, double max_dist, double& dist);

        template<class T2>
        bool    test_equal(const T2& res, const T2& res_gen, double max_dist, double& dist);

        template<class T, class TMP>
        bool    test_equal_fma(int size, const T* res, const TMP* res_gen, double max_dist, double& dist);

        template<class T, class TMP>
        bool    test_equal_fma(const T& res, const TMP& res_gen, double max_dist, double& dist);

        template<class Func, class T2>
        bool    test_constraints_1(int size, const T2* res, const T2* in);

        template<class Func, class T2>
        bool    test_constraints_2(int size, const T2* res, const T2* in1, const T2* in2);

        template<class Func, class T2>
        bool    test_constraints_elem_1(const T2& res, const T2& in1);

        template<class Func, class T2>
        bool    test_constraints_elem_2(const T2& res, const T2& in1, const T2& in2);

        template<class T2>
        bool    is_normalized(const T2& res);

        template<class T2, class TMP, class Func>
        double  test_function_gen(int size, int n_rep, const T2* in, TMP* out);

        template<class T, class T2, class Func>
        double  test_function_base(int size, int n_rep, const T2* in, T2* out);

        template<class T2, class Func>
        double  test_function_twofold(int size, int n_rep, const T2* in, T2* out);

        template<class T2, class TMP, class Func>
        double  test_function_bin_gen(int size, int n_rep, const T2* in1, const T2* in2, TMP* out);

        template<class T, class T2, class Func>
        double  test_function_bin_base(int size, int n_rep, const T2* in1, const T2* in2, T2* out);

        template<class T2, class Func>
        double  test_function_bin_twofold(int size, int n_rep, const T2* in1, const T2* in2, T2* out);

        template<class T, class TMP, class Func>
        double  test_function_fma_gen(int size, int n_rep, const T* in1, const T* in2, 
                    const T* in3, TMP* out);

        template<class T, class Func>
        double  test_function_fma_base(int size, int n_rep, const T* in1, const T* in2, 
                    const T* in3, T* out);

        template<class T, class Func>
        double  test_function_fma_twofold(int size, int n_rep, const T* in1, const T* in2, 
                    const T* in3, T* out);
};

}}
