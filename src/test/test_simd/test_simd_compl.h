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

#include <string>

namespace matcl { namespace test
{

void test_performance_complex();
void test_values_complex();

class test_simd_compl
{
    private:
        std::string m_instr_tag;
        bool        m_test_values;

    public:
        test_simd_compl(bool test_values);

        void    make_binary();
        void    make_unary();

    private:
        int     get_size() const;
        int     get_size_perf() const;
        int     get_num_rep() const;

        template<class T>
        bool    test_equal(int size, const T* res, const T* res_gen, double max_dist, double& dist,
                    bool componentwise);

        template<class T>
        bool    test_equal(const T& res, const T& res_gen, double max_dist, double& dist,
                    bool componentwise);

        template<class T>
        bool    test_equal_real(const T& res, const T& res_gen, double max_dist, double& dist);

        template<class T>
        void    test_functions();

        template<class T>
        void    test_functions_bin();

        template<class T, class Func>
        void    test_function(formatted_disp& fd, int size, const T* in, 
                    T* out, T* out_gen, bool test_componentwise);

        template<class T, class Func>
        void    test_function_block(formatted_disp& fd, int size, const T* in, 
                    T* out, T* out_gen, bool test_componentwise);

        template<class T, class Func>
        void    test_function_bin(formatted_disp& fd, int size, const T* in_1, 
                    const T* in_2, T* out, T* out_gen, double max_dist, bool test_componentwise);

        template<class T, class Func, class TR = typename ms::details::real_type<T>::type>
        void    test_function_bin_RC(formatted_disp& fd, int size, const TR* in_1, 
                    const T* in_2, T* out, T* out_gen, double max_dist, bool test_componentwise);

        template<class T, class Func, class TR = typename ms::details::real_type<T>::type>
        void    test_function_bin_CR(formatted_disp& fd, int size, const T* in_1, 
                    const TR* in_2, T* out, T* out_gen, double max_dist, bool test_componentwise);

        template<class T, class Simd_type, class Func>
        double  test_function_simd(int size, int n_rep, const T* in, T* out);

        template<class T, class Simd_type, class Func>
        double  test_function_bin_simd(int size, int n_rep, const T* in1, const T* in2, T* out);

        template<class T, class TR, class Simd_type, class Func>
        double  test_function_bin_simd_RC(int size, int n_rep, const TR* in1, const T* in2, T* out);

        template<class T, class TR, class Simd_type, class Func>
        double  test_function_bin_simd_CR(int size, int n_rep, const T* in1, const TR* in2, T* out);

        std::string make_status(bool ok, double ulp_dist);
};

}}
