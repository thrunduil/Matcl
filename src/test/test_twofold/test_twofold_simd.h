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

void test_functions_simd();
void test_io_simd();
void test_error_simd();

class test_twofold_simd
{
    public:
        void    make_unary();
        void    make_binary();
        void    make_io();
        void    make_error();

    private:
        int     get_size() const;
        int     get_size_perf() const;
        int     get_num_rep() const;

        template<class T>
        std::string get_header() const;

        template<class T>
        void    test_functions();

        template<class T>
        void    test_functions_bin();

        template<class T>
        void    test_functions_io();

        template<class Float_type>
        void    make_error_type();

        template<class Simd>
        bool    make_error(const twofold<Simd>& x, int mult);

        template<class T, class Func>
        void    test_function(formatted_disp& fd, int size, const T* in_1, 
                    const T* in_2, T* out_1, T* out_2, T* out_1_base);

        template<class T, class Func>
        void    test_function_io(formatted_disp& fd, int size, const T* in_1, 
                    const T* in_2, T* out_1, T* out_2, T* out_1_base);

        template<class T, class TS, class Func>
        double  test_function_base(int size, int n_rep, const T* in1, T* out_1);

        template<class T, class TS, class Func>
        double  test_function_twofold(int size, int n_rep, const T* in_1, const T* in_2, 
                    T* out_1, T* out_2);

        template<class T>
        void    test_functions_bin(formatted_disp& fd, int size, T* in_11, 
                    T* in_12, T* in_21, T* in_22, T* out_1, T* out_2, T* out_base);

        template<class T>
        void    test_functions_sum_bin(formatted_disp& fd, int size, T* in_11, 
                    T* in_12, T* in_21, T* in_22, T* out_1, T* out_2, T* out_base);

        template<class T, class Func>
        void    test_function_bin(formatted_disp& fd, int size, const T* in_11, 
                    const T* in_12, const T* in_21, const T* in_22, T* out_1, T* out_2, 
                    T* out_base, bool inaccurate);

        template<class T, class TS, class Func>
        double  test_function_bin_base(int size, int n_rep, const T* in1, const T* in2, T* out);

        template<class T, class TS, class Func>
        double  test_function_bin_twofold(int size, int n_rep, const T* in_11, const T* in_12, 
                    const T* in_21, const T* in_22, T* out_1, T* out_2);

        template<class T>
        bool    test_equal(int size, const T* res, const T* res_gen, double max_dist, double& dist);

        template<class T>
        bool    test_equal_inacc(int size, const T* in_11, const T* in_21, const T* res, 
                    const T* res_gen, double max_dist, double& dist);

        template<class T>
        bool    test_equal(const T& res, const T& res_gen, double max_dist, double& dist);

        template<class T>
        bool    test_equal_inacc(const T& in_11, const T& in_21, const T& res, const T& res_gen, 
                    double max_dist, double& dist);
};

}}
